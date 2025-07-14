/*!
  \file partition_view.hpp
  \brief Implements partition interface for AIG network

  \author Jingren Wang
*/

#include "fanout_view.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <mockturtle/networks/aig.hpp>

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

  explicit partition_view( aig_network const& ntk, partition_view_params const& ps = {} ) : _hypNum( 0 ), _ps( ps )
  {
    collect_hypgraph( ntk );
    write_hypgraph( ntk, ps );
  }

  ~partition_view()
  {
  }

private:
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
};

} // namespace mockturtle
