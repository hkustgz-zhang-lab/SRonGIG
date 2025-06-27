#include "experiments.hpp"
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>

int main()
{
  /*
  Note:
  In MIG, create_and and create_not method are also implemented, so here directly reading an AIG file to MIG get same node number and same level number. But actually MIG node number can be smaller, check how experiments/mapper.cpp does it.
  MIG converts back is difference, node number may increase, since the process is convert to KLUT network first and then back to aig, there's possiblility of loss.
  */

  using namespace experiments;
  using namespace mockturtle;
  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, bool> exp( "readerLoss", "benchmark", "nodes num before", "nodes num after", "num back to aig", "depth before", "depth after", "level back to aig", "cec res" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    mig_network mig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t aig_size_before = aig.num_gates();
    const uint32_t aig_depth_before = depth_view( aig ).depth();

    const uint32_t mig_size_after = mig.num_gates();
    const uint32_t mig_depth_after = depth_view( mig ).depth();

    klut_network klut = lut_map( mig );
    aig_network aig_back = convert_klut_to_graph<aig_network>( klut );

    const uint32_t aig_back_size = aig_back.num_gates();
    const uint32_t aig_back_depth = depth_view( aig_back ).depth();

    const auto cec1 = benchmark == "hyp" ? true : abc_cec( mig, benchmark );

    exp( benchmark, aig_size_before, mig_size_after, aig_back_size, aig_depth_before, mig_depth_after, aig_back_depth, cec1 );
  }
  exp.save();
  exp.table();
  return 0;
}