#include <catch.hpp>
#include <fmt/core.h>
#include <mtkahypar.h>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/views/partition_view.hpp>

using namespace mockturtle;

TEST_CASE( "create and dump hypergraph from an AIG", "[partition]" )
{
  aig_network aig;
  const auto x1 = aig.create_pi();
  const auto x2 = aig.create_pi();
  const auto x3 = aig.create_pi();
  const auto x4 = aig.create_pi();

  const auto f1 = aig.create_and( x1, x2 );
  const auto f2 = aig.create_and( x3, x4 );
  const auto f3 = aig.create_and( x1, x3 );
  const auto f4 = aig.create_and( f1, f2 );
  const auto f5 = aig.create_and( f3, f4 );

  aig.create_po( f5 );

  partition_view_params ps;
  ps.file_name = fmt::format( "{}/test.hmetis", PARTITION_TEST_PATH );
  partition_view aig_p{ aig, ps };

  ps.si_w_on_hyperedges = true;
  ps.file_name = fmt::format( "{}/test_edge_weight.hmetis", PARTITION_TEST_PATH );
  partition_view aig_p_e_w{ aig, ps };

  ps.si_w_on_hyperedges = false;
  ps.si_w_on_vertices = true;
  ps.file_name = fmt::format( "{}/test_vertices_weight.hmetis", PARTITION_TEST_PATH );
  partition_view aig_p_v_w{ aig, ps };

  ps.si_w_on_hyperedges = true;
  ps.si_w_on_vertices = true;
  ps.file_name = fmt::format( "{}/test_edge_vertices_weight.hmetis", PARTITION_TEST_PATH );
  partition_view aig_p_e_v_w{ aig, ps };
}
