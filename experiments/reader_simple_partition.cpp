#include "experiments.hpp"
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mtkahypar.h>
#include <thread>

int main()
{
  /*
  Note:
  In MIG, create_and and create_not method are also implemented, so here directly reading an AIG file to MIG get same node number and same level number. But actually MIG node number can be smaller, check how experiments/mapper.cpp does it.
  MIG converts back is difference, node number may increase, since the process is convert to KLUT network first and then back to aig, there's possiblility of loss.
  */

  using namespace experiments;
  using namespace mockturtle;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    // TODO: @Jingren Read from mt directly
    // aig_network aig;
    // if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    //{
    //   continue;
    // }

    // const uint32_t aig_size_before = aig.num_gates();
    // const uint32_t aig_depth_before = depth_view( aig ).depth();

    std::string command = fmt::format( "abc -q \"read_aiger {}; write_hmetis /tmp/adder.hmetis\"", benchmark_path( "adder" ) );
#if WIN32
    std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
    std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
#endif

    // use the official repo's example as a test
    mt_kahypar_error_t error{};

    // Initialize
    mt_kahypar_initialize(
        std::thread::hardware_concurrency() /* use all available cores */,
        true /* activate interleaved NUMA allocation policy */ );

    // Setup partitioning context
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset( DEFAULT );
    // In the following, we partition a hypergraph into two blocks
    // with an allowed imbalance of 3% and optimize the connective metric (KM1)
    mt_kahypar_set_partitioning_parameters( context,
                                            2 /* number of blocks */, 0.03 /* imbalance parameter */,
                                            KM1 /* objective function */ );
    mt_kahypar_set_seed( 42 /* seed */ );
    // Enable logging
    mt_kahypar_status_t status =
        mt_kahypar_set_context_parameter( context, VERBOSE, "1", &error );
    assert( status == SUCCESS );

    // Load Hypergraph for DEFAULT preset
    mt_kahypar_hypergraph_t hypergraph =
        mt_kahypar_read_hypergraph_from_file( "/tmp/adder.hmetis",
                                              context, HMETIS /* file format */, &error );
    if ( hypergraph.hypergraph == nullptr )
    {
      std::cout << error.msg << std::endl;
      std::exit( 1 );
    }

    // Partition Hypergraph
    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
        mt_kahypar_partition( hypergraph, context, &error );
    if ( partitioned_hg.partitioned_hg == nullptr )
    {
      std::cout << error.msg << std::endl;
      std::exit( 1 );
    }

    // Extract Partition
    auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(
        mt_kahypar_num_hypernodes( hypergraph ) );
    mt_kahypar_get_partition( partitioned_hg, partition.get() );

    // Extract Block Weights
    auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>( 2 );
    mt_kahypar_get_block_weights( partitioned_hg, block_weights.get() );

    // Compute Metrics
    const double imbalance = mt_kahypar_imbalance( partitioned_hg, context );
    const int km1 = mt_kahypar_km1( partitioned_hg );

    // Output Results
    std::cout << "Partitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;

    mt_kahypar_free_context( context );
    mt_kahypar_free_hypergraph( hypergraph );
    mt_kahypar_free_partitioned_hypergraph( partitioned_hg );
  }
  return 0;
}