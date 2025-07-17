/*!
  \file partition_view.hpp
  \brief Implements partition interface for AIG network

  \author Jingren Wang
*/

#include "fanout_view.hpp"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/network_utils.hpp>
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

  /*! \brief Write out to the specific file name */
  std::string file_name{ "tmp.hmetis" };
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

  explicit partition_view( aig_network const& ntk, partition_view_params const& ps = {} ) : _hypNum( 0 ), _ps( ps ), _ntk( ntk )
  {
    collect_hypgraph( ntk );
    write_hypgraph( ntk, ps );
  }

  ~partition_view()
  {
  }

  std::vector<std::tuple<aig_network, std::vector<node>, std::vector<signal>>> construct_from_partition( int nPart, std::map<edge_id, std::vector<node_id_k>> vBoundaries, std::map<node_id_k, block_id> const& node_block_io, std::map<node_id_k, block_id>& node_block )
  {
    // create nPart aig_network in parallel
    std::vector<std::tuple<aig_network, std::vector<node>, std::vector<signal>>> vAigs_win( nPart );
    // construct aig from scratch
    std::cout << "Size of the io map is " << node_block_io.size() << std::endl;
    /*
    block0: vCis; vCos; vNodes;
    block1: vCis; vCos; vNodes;
    ...
    */
    std::vector<std::vector<node>> node_id_split( 3 * nPart );

    /*
    For all pis/pos, these can only be split, even after patition, they are still pis/pos ,don't need to check the boundaries, just might belongs to different block.
    */
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
        po_index = _ntk.num_gates() + _ntk.num_pis();
      }
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
    } );

    /*
    For each boundary, we should identify the source and sinks:
    There's always a SINGLE source and maybe multiple sinks, this is by design of the AIG, so for a group of nodes attached to a single hyperEdge, the one has the lowest level is the SOURCE.
    */
    depth_view<aig_network> ntk_dep{ _ntk };
    for ( auto& edg : vBoundaries )
    {
      assert( edg.second.size() != 0 );
      auto src_index = edg.second[0];
      if ( src_index == ( _ntk.num_gates() + _ntk.num_pis() ) )
      {
        src_index = 0;
      }
      auto node_source = _ntk.index_to_node( src_index );
      auto index = 1;
      while ( index < edg.second.size() )
      {
        auto target_to_compare = edg.second[index];
        if ( target_to_compare == _ntk.num_gates() + _ntk.num_pis() )
        {
          assert( target_to_compare == _ntk.num_gates() + _ntk.num_pis() );
          std::cout << "Encounter with const0 id " << target_to_compare << std::endl;
          target_to_compare = 0;
        }
        if ( ntk_dep.level( _ntk.index_to_node( target_to_compare ) ) < ntk_dep.level( node_source ) )
        {
          node_source = _ntk.index_to_node( target_to_compare );
        }
        index++;
      }
      // Add SOURCE to the fanout of the block itself belongs to
      auto b_s = _ntk.node_to_index( node_source );
      if ( b_s == 0 )
      {
        b_s = _ntk.num_gates() + _ntk.num_pis();
      }
      auto& block_id_s = node_block[b_s];
      auto& ac_pos = node_id_split[block_id_s * 3 + 1];
      // only single output is needed(?)
      if ( !check_node_exist( node_id_split[block_id_s * 3 + 1], node_source ) )
      {
        ac_pos.push_back( node_source );
      }

      // SOURCE now is node_source
      /*
      Use fanout iteration method, check:
      The number of different id block in fanouts of the SOURCE show how many CI should be provided for the according block.
      */
      index = 0;
      while ( index < edg.second.size() )
      {
        auto it = node_block.find( edg.second[index] );
        // if exist and not the source of the hyperEdge
        if ( it != node_block.end() && edg.second[index] != _ntk.node_to_index( node_source ) )
        {
          // if node in fanout has different block id as SOURCE, then add SOURCE to CI of the fanout
          auto tmp_fo_block = node_block[edg.second[index]];
          if ( tmp_fo_block != node_block[_ntk.node_to_index( node_source )] )
          {
            // check if the block already conatins the CI as node_source
            auto& ac_pis = node_id_split[tmp_fo_block * 3];
            if ( check_node_exist( ac_pis, _ntk.index_to_node( node_source ) ) )
            {
              index++;
              continue;
            }
            else
            {
              ac_pis.push_back( _ntk.index_to_node( node_source ) );
            }
          }
          index++;
        }
        else
        {
          index++;
          continue;
        }
      }
    }

    auto ori_num_gate = _ntk.num_gates();
    std::cout << "Num of gate in original ntk : " << ori_num_gate << std::endl;
    _ntk.foreach_gate( [&]( auto const& n_and ) {
      auto and_index = _ntk.node_to_index( n_and );
      auto block_id_and = node_block[and_index];
      auto cis_id = node_id_split[block_id_and * 3];
      // Do not contain PIs, as suggested by clone_subnetwork api
      if ( !check_node_exist( cis_id, n_and ) )
      {
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
      std::cout << "================\nBlock " << i << std::endl;
      aig_network win;
      std::vector<node> inputs = node_id_split[i * 3];
      assert( inputs.size() > 0 );
      std::cout << "Size of CI: " << inputs.size() << std::endl;
      std::vector<node> outputs_o = node_id_split[i * 3 + 1];
      assert( outputs_o.size() > 0 );
      std::vector<signal> outputs;
      construct_out_sig( outputs_o, outputs );
      std::cout << "Size of CO: " << outputs.size() << std::endl;
      std::vector<node> gates = node_id_split[i * 3 + 2];
      std::cout << "Size of gates: " << gates.size() << std::endl;
      count_all_gate += gates.size();
      clone_subnetwork( _ntk, inputs, outputs, gates, win );
      std::get<0>( vAigs_win[i] ) = win;
      std::get<1>( vAigs_win[i] ) = inputs;
      std::get<2>( vAigs_win[i] ) = outputs;
    }
    assert( count_all_gate == ori_num_gate );
    return vAigs_win;
  }

private:
  void construct_out_sig( std::vector<node> const& outputs_o, std::vector<signal>& outputs )
  {

    for ( auto& outs_o : outputs_o )
    {
      _ntk.foreach_node( [&]( auto const& n ) {
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          if ( _ntk.get_node( f ) == outs_o )
          {
            outputs.push_back( _ntk.make_signal( outs_o ) );
          }
        } );
      } );
    }
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
    _sinkHyp.erase( 0 );
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
    os << _hypNum << " " << ntk.num_gates() + ntk.num_pis() << "\n";
    for ( auto item = _sinkHyp.begin(); item != _sinkHyp.end(); item++ )
    {

      if ( item->second.size() == 0 )
      {
        continue;
      }
      else
      {
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
    os << "%% Mockturtle finished writing the hMetis file." << std::endl;
    os.close();
  }
  int _hypNum;
  partition_view_params _ps;
  std::map<node, std::vector<node>> _sinkHyp;
  aig_network _ntk;
};

} // namespace mockturtle
