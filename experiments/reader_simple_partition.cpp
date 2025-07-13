#include "experiments.hpp"
#include <_stdio.h>
#include <fmt/core.h>
#include <lorina/aiger.hpp>
#include <memory>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mtkahypar.h>
#include <thread>

#include <iostream>   // For std::cout and std::cerr
#include <sys/stat.h> // For mkdir

#include <filesystem> // Requires C++17
#include <iostream>

int creatDir( const std::filesystem::path* dir_path )
{

  std::error_code ec; // For error handling
  if ( std::filesystem::exists( *dir_path ) )
  {
    std::cout << "Directory exists." << std::endl;
    return 0;
  }
  if ( std::filesystem::create_directory( *dir_path, ec ) )
  {
    std::cout << "Directory created successfully: " << *dir_path << std::endl;
    return 0;
  }
  else
  {
    std::cerr << "Failed to create directory: " << *dir_path << std::endl;
    std::cerr << "Error: " << ec.message() << std::endl;
    return 1;
  }
  return 0;
}

int main()
{
  /*
  Note:
  Due to lack of proper documentation, using mt-kapypar is time-consuming, the following cases give some simple examples.
  */
  std::cout << EXPERIMENTS_PATH << std::endl;
  std::filesystem::path pHyOut = fmt::format( "{}HypOut", EXPERIMENTS_PATH );
  const std::string pHyOutS = pHyOut;

  if ( creatDir( &pHyOut ) == 0 )
  {
    std::cout << "Successfully created the output directory " << pHyOut << std::endl;
  }
  else
  {
    std::exit( 1 );
  }

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

    std::string command = fmt::format( "abc -q \"read_aiger {}; write_hmetis {}/tmp.hmetis;\"", benchmark_path( benchmark ), pHyOutS );
#if WIN32
    std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
#else
    // std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
#endif
    // Use system call here, the above method leads to IO error
    system( command.c_str() );

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
                                            3 /* number of blocks */, 0.03 /* imbalance parameter */,
                                            KM1 /* objective function */ );
    mt_kahypar_set_seed( 42 /* seed */ );
    // Enable logging
    mt_kahypar_status_t status =
        mt_kahypar_set_context_parameter( context, VERBOSE, "0", &error );
    assert( status == SUCCESS );

    // Load Hypergraph for DEFAULT preset
    mt_kahypar_hypergraph_t hypergraph =
        mt_kahypar_read_hypergraph_from_file( fmt::format( "{}/tmp.hmetis", pHyOutS ).c_str(),
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
    auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>( 3 );
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
    std::cout << "Weight of Block 2 = " << block_weights[1] << std::endl;
    std::cout << std::endl;

    // You cannot give a directory that doesn't exist here, std::ofstream doesn't create new directory by default
    mt_kahypar_status_t status_out = mt_kahypar_write_partition_to_file( partitioned_hg,
                                                                         fmt::format( "{}/{}hypOut", pHyOutS, benchmark ).c_str(),
                                                                         &error );
    assert( status_out == SUCCESS );
    if ( remove( fmt::format( "{}/tmp.hmetis", pHyOutS ).c_str() ) == 0 )
    {
      std::cout << "Successfully deleted the tmp file." << std::endl;
    }
    else
    {
      std::cout << "Failed to delete the tmp file." << std::endl;
    }

    mt_kahypar_free_context( context );
    mt_kahypar_free_hypergraph( hypergraph );
    mt_kahypar_free_partitioned_hypergraph( partitioned_hg );
  }
  return 0;
}