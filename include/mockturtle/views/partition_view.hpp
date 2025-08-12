/*!
  \file partition_view.hpp
  \brief Implement partition interface for AIG network

  \author Jingren Wang
*/

#include "fanout_view.hpp"
#include "mtkahypar.h"
#include "mtkahypartypes.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/utils/network_utils.hpp>
#include <mockturtle/utils/window_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <ostream>
#include <string>
#include <vector>

namespace mockturtle
{

struct partition_view_params
{
  /*! \brief Skip po as sink explictly. */
  bool skip_po_as_sink{ true };

  /*! \brief Simulate weight on edges. */
  bool si_w_on_hyperedges{ false };

  /*! \brief Simulate weight on vertices. */
  bool si_w_on_vertices{ false };

  /*! \brief Write out to the specific file name */
  std::string file_name{ "tmp.hmetis" };

  /*! \brief Number of partition block. */
  mt_kahypar_partition_id_t num_blocks{ 2 };

  /*! \brief Seed of the function. */
  size_t seed{ 42 };

  /*! \brief Epsilon */
  double epsilon{ 0.03 };
};

class partition_view
{
public:
  using storage = typename aig_network::storage;
  using node = typename aig_network::node;
  using signal = typename aig_network::signal;
  using node_id_k = unsigned long int;
  using block_id = int;
  using edge_id = unsigned long int;

  explicit partition_view( aig_network const& ntk, partition_view_params const& ps = {} ) : _hypNum( 0 ), _ps( ps ), _ntk( ntk ), refs( ntk.size() )
  {
    _sinkHyp.clear();
    collect_hypgraph( ntk );
    write_hypgraph( ntk, ps );
  }

  ~partition_view()
  {
  }

  std::vector<std::tuple<aig_network, std::vector<node>, std::vector<signal>, std::vector<node>>> construct_from_partition( int nPart, const std::unique_ptr<mt_kahypar_partition_id_t[]>& partition, const mt_kahypar_hypergraph_t& hypergraph )
  {
    auto node_block = convertToMap( partition, mt_kahypar_num_hypernodes( hypergraph ) );
    // create nPart aig_network in parallel
    std::vector<std::tuple<aig_network, std::vector<node>, std::vector<signal>, std::vector<node>>> vAigs_win( nPart );
    // construct aig from scratch
    /*
    block0: vCis; vCos; vNodes;
    block1: vCis; vCos; vNodes;
    ...
    */
    std::vector<std::vector<node>> node_id_split( 3 * nPart );

    /*
    For all pis/pos, these can only be split, even after patition, they are still pis/pos, don't need to check the boundaries, just might belongs to different block.
    */
    // Deal with constant0 first
    /*
    Eg Inputs 1,2,3 ANDs 4,5,6, then _ntk.num_gates() + _ntk.num_pis() + 1 = 7
    If _ntk.num_gates() + _ntk.num_pis() + 1 means no zero in hMetis
    Eg Inputs 1,2,3 ANDs 4,5,6, Zero 7, then _ntk.num_gates() + _ntk.num_pis() + 2 = 8
    _ntk.num_gates() + _ntk.num_pis() + 2 means zero is counted.

    Currently not supporting const zero.
    */
    assert( node_block.size() == _ntk.num_gates() + _ntk.num_pis() + 1 );
    // auto block_id_const = node_block[_ntk.num_gates() + _ntk.num_pis() + 1];
    // node_id_split[block_id_const * 3].push_back( _ntk.get_node(_ntk.get_constant(false)) );

    _ntk.foreach_pi( [&]( auto const& n_pi ) {
      auto pi_index = _ntk.node_to_index( n_pi );
      auto it = node_block.find( pi_index );
      if ( it == node_block.end() )
      {
        std::cout << "Something wrong. There's no PI index " << pi_index << std::endl;
      }
      else
      {
        auto block_id = it->second;
        node_id_split[block_id * 3].push_back( n_pi );
      }
    } );
    _ntk.foreach_po( [&]( auto const& n_po ) {
      auto po_index = _ntk.node_to_index( _ntk.get_node( n_po ) );
      // should be one of the and type or PI type
      assert( _ntk.is_and( _ntk.get_node( n_po ) ) || _ntk.is_pi( _ntk.get_node( n_po ) ) || _ntk.is_constant( _ntk.get_node( n_po ) ) );
      if ( po_index == 0 )
      {
        std::cout << "[Warn] PO has const 0." << std::endl;
        /*
        @TODO Jingren Wang
        */
      }
      else
      {
        auto it = node_block.find( po_index );
        if ( it == node_block.end() )
        {
          std::cout << "Something wrong. There's no PO index " << po_index << std::endl;
        }
        else
        {
          auto block_id = it->second;
          node_id_split[block_id * 3 + 1].push_back( _ntk.get_node( n_po ) );
        }
      }
    } );

    auto ori_num_gate = _ntk.num_gates();
    std::cout << "Num of gate in original ntk : " << ori_num_gate << std::endl;
    _ntk.foreach_gate( [&]( auto const& n_and ) {
      auto and_index = _ntk.node_to_index( n_and );
      auto block_id_and = node_block[and_index];
      auto cis_id = node_id_split[block_id_and * 3];
      // Do not contain PIs, as suggested by clone_subnetwork api
      if ( !check_node_exist( cis_id, n_and ) )
      {
        assert( _ntk.node_to_index( n_and ) <= _ntk.num_pis() + _ntk.num_gates() );
        node_id_split[block_id_and * 3 + 2].push_back( n_and );
      }
    } );
    std::cout << "Writing out aigs..." << std::endl;
    /*
    Use window based method to construct a subnetwork and each one of them can be insert back to original aig network and maintain equivalence
    */
    auto count_all_gate = 0;
    for ( auto i = 0; i < nPart; i++ )
    {
      aig_network win;
      color_view c_ntk{ _ntk };
      std::vector<node> gates = node_id_split[i * 3 + 2];
      std::vector<node> inputs_w = collect_inputs( c_ntk, gates );
      assert( inputs_w.size() > 0 );
      std::stable_sort( std::begin( inputs_w ), std::end( inputs_w ) );
      std::stable_sort( std::begin( gates ), std::end( gates ) );
      std::vector<signal> outputs_w = collect_outputs( c_ntk, inputs_w, gates, refs );
      assert( outputs_w.size() > 0 );
      count_all_gate += gates.size();
      clone_subnetwork( _ntk, inputs_w, outputs_w, gates, win );
      std::get<0>( vAigs_win[i] ) = win;
      std::get<1>( vAigs_win[i] ) = inputs_w;
      std::get<2>( vAigs_win[i] ) = outputs_w;
      std::get<3>( vAigs_win[i] ) = gates;
    }
    assert( count_all_gate == ori_num_gate );
    return vAigs_win;
  }

  void insert_back( std::tuple<aig_network, std::vector<node>, std::vector<signal>, std::vector<node>> const& aig_part )
  {
    auto aig = _ntk;
    // create signals
    std::vector<aig_network::signal> i_sigs;
    for ( auto const& i : std::get<1>( aig_part ) )
    {
      i_sigs.push_back( aig.make_signal( i ) );
    }
    uint32_t counter = 0u;
    // insert back now
    color_view c_aig{ aig };
    assert( count_reachable_dead_nodes( c_aig ) == 0u );
    insert_ntk( aig, i_sigs.begin(), i_sigs.end(), std::get<0>( aig_part ), [&]( aig_network::signal const& _new ) {
      assert( !c_aig.is_dead( c_aig.get_node( _new ) ) );
      auto const _old = std::get<2>( aig_part ).at( counter++ );
      if ( _old == _new )
      {
        return true;
      }

      if ( _old != _new )
      {
        aig.substitute_node( aig.get_node( _old ), aig.is_complemented( _old ) ? !_new : _new );
      }
      return true;
    } );
  }

private:
  std::map<mt_kahypar_hypernode_id_t, mt_kahypar_partition_id_t> convertToMap(
      const std::unique_ptr<mt_kahypar_partition_id_t[]>& array,
      size_t size )
  {
    std::map<mt_kahypar_hypernode_id_t, mt_kahypar_partition_id_t> result;
    result[0] = -1;
    for ( size_t i = 0; i < size; ++i )
    {
      result[i + 1] = array[i];
    }
    assert( result.size() == size + 1 );
    return result;
  }

  bool check_node_exist( std::vector<node> const& nodes, node const& p )
  {
    for ( auto& n : nodes )
    {
      if ( n == p )
      {
        return true;
      }
    }
    return false;
  }

  void collect_single_node( fanout_view<aig_network> const& ntk, node const& n )
  {
    _sinkHyp[n] = {};
    ntk.foreach_fanout( n, [&]( auto const& fo ) {
      _sinkHyp[n].push_back( fo );
    } );
  }

  void collect_hypgraph( aig_network const& ntk )
  {
    fanout_view f_aig{ ntk };
    // don't need to consider this since mt doesn't create new node for PO
    if ( !_ps.skip_po_as_sink )
    {
    }
    ntk.foreach_pi( [&]( auto pi ) {
      collect_single_node( f_aig, pi );
    } );
    ntk.foreach_node( [&]( auto node ) {
      collect_single_node( f_aig, node );
    } );
    std::vector<node> copy = _sinkHyp[0];
    auto last_index = _sinkHyp.size();
    _sinkHyp[0] = {};
    _sinkHyp[last_index] = copy;
  }

  void hypNumCount( std::map<node, std::vector<node>> const& allHyp, int& hypNum )
  {
    for ( auto& item : allHyp )
    {
      if ( item.second.size() > 0 )
      {
        hypNum++;
      }
    }
  }

  void write_hypgraph( aig_network const& ntk, partition_view_params ps )
  {
    // skip the weight simulation for now
    std::ofstream os( ps.file_name.c_str(), std::ofstream::out );
    hypNumCount( _sinkHyp, _hypNum );
    bool constExist = false;
    if ( _sinkHyp[_sinkHyp.size() - 1].size() == 0 )
    {
      os << _hypNum << " " << ntk.num_gates() + ntk.num_pis();
    }
    else
    {
      os << _hypNum << " " << ntk.num_gates() + ntk.num_pis() + 1;
      constExist = true;
    }

    if ( ps.si_w_on_hyperedges && ps.si_w_on_vertices )
    {
      os << " " << 11;
    }
    else if ( ps.si_w_on_hyperedges )
    {
      os << " " << 1;
    }
    else if ( ps.si_w_on_vertices )
    {
      os << " " << 10;
    }
    os << "\n";

    for ( auto item = _sinkHyp.begin(); item != _sinkHyp.end(); item++ )
    {

      if ( item->second.size() == 0 )
      {
        continue;
      }
      else
      {
        if ( ps.si_w_on_hyperedges )
        {
          // This could be customized, currently use the fanout number as the weight
          os << item->second.size() << " ";
        }
        os << ntk.node_to_index( item->first ) << " ";
        for ( auto t = item->second.begin(); t < item->second.end(); t++ )
        {
          if ( t == item->second.end() - 1 )
          {
            os << ntk.node_to_index( *t ) << "\n";
          }
          else
          {
            os << ntk.node_to_index( *t ) << " ";
          }
        }
      }
    }
    // simulate weights of vertices to 1, which could be modified
    if ( ps.si_w_on_vertices )
    {
      for ( auto item = _sinkHyp.begin(); item != _sinkHyp.end(); item++ )
      {
        if ( item->second.size() == 0 )
        {
          continue;
        }
        os << 1 << "\n";
      }
    }
    os << "%% Mockturtle finished writing the hMetis file." << std::endl;
    if ( constExist )
    {
      // This should not be triggered currently without const zero.
      assert( 0 );
      os << "%% Const exists as the largest index." << std::endl;
    }
    os.close();
  }
  int _hypNum;
  partition_view_params _ps;
  std::map<node, std::vector<node>> _sinkHyp;
  aig_network _ntk;
  std::vector<uint32_t> refs;
};

} // namespace mockturtle
