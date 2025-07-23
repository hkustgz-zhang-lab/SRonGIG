#include "experiments.hpp"
#include <_stdio.h>
#include <cassert>
#include <cstddef>
#include <fmt/core.h>
#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <memory>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/partition_view.hpp>
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
    if ( benchmark == "hyp" )
    {
      continue;
    }
    fmt::print( "[i] processing {}\n", benchmark );

    // TODO: @Jingren Read from mt directly
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    auto ori_gate_num = aig.num_gates();

    partition_view_params ps;
    ps.file_name = fmt::format( "{}/tmp_ctrl.hmetis", pHyOutS );
    partition_view aig_p{ aig, ps };

    // Althernative ABC method
    //    std::string command = fmt::format( "abc -q \"read_aiger {}; write_hmetis {}/tmp.hmetis;\"", benchmark_path( benchmark ), pHyOutS );
    // #if WIN32
    //    std::unique_ptr<FILE, decltype( &_pclose )> pipe( _popen( command.c_str(), "r" ), _pclose );
    // #else
    //    // std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
    // #endif
    //    // Use system call here, the above method leads to IO error
    //    system( command.c_str() );

    // use the official repo's example as a test
    mt_kahypar_error_t error{};

    // Initialize
    mt_kahypar_initialize(
        std::thread::hardware_concurrency() /* use all available cores */,
        true /* activate interleaved NUMA allocation policy */ );

    // Setup partitioning context
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset( DETERMINISTIC );
    // In the following, we partition a hypergraph into two blocks
    // with an allowed imbalance of 3% and optimize the connective metric (KM1)
    mt_kahypar_set_partitioning_parameters( context,
                                            2 /* number of blocks */, 0.03 /* imbalance parameter */,
                                            KM1 /* objective function */ );
    mt_kahypar_set_seed( 42 /* seed */ );
    // Enable logging
    mt_kahypar_status_t status =
        mt_kahypar_set_context_parameter( context, VERBOSE, "0", &error );
    assert( status == SUCCESS );

    // Load Hypergraph for DEFAULT preset
    mt_kahypar_hypergraph_t hypergraph =
        mt_kahypar_read_hypergraph_from_file( fmt::format( "{}/tmp_ctrl.hmetis", pHyOutS ).c_str(),
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

    // Iterate over all edges
    /*
    The idea is that the connectivity larger than 1 is identified as crossing multiple blocks.
    Get the boundaries and those nodes connected with boundaries.
    There's also a problem about the direction, let's hope after the depth view this can be solved easily.
    */

    std::map<mt_kahypar_hyperedge_id_t, std::vector<mt_kahypar_hypernode_id_t>> vBoundaries;
    for ( mt_kahypar_hyperedge_id_t edge = 0; edge < mt_kahypar_num_hyperedges( hypergraph ); ++edge )
    {
      // std::cout << "Edge " << edge << std::endl;
      //  Print number of blocks connected by edge
      mt_kahypar_partition_id_t connectedBlockNum = mt_kahypar_connectivity( partitioned_hg, edge );
      // std::cout << "Connectivity = " << connectedBlockNum << std::endl;
      std::vector<mt_kahypar_hypernode_id_t> pins_buffer;
      pins_buffer.resize( mt_kahypar_hyperedge_size( hypergraph, edge ) );
      mt_kahypar_get_hyperedge_pins( hypergraph, edge, pins_buffer.data() );
      // mt_kahypar_get_hyperedge_pins(hypergraph, edge, pins_buffer.get());
      if ( connectedBlockNum > 1 )
      {
        if ( vBoundaries.count( edge ) > 0 )
        {
          std::cout << "Encounter with edges twice, this shouldn't happen." << std::endl;
        }
        else
        {
          vBoundaries[edge] = pins_buffer;
        }
      }
    }
    std::cout << "Number of the boundaries " << vBoundaries.size() << std::endl;
    // show the map structure in detail
    for ( auto& item : vBoundaries )
    {
      mt_kahypar_hypernode_id_t hypId = 0;
      while ( hypId < item.second.size() )
      {
        // Recover the original id in hMetis
        item.second[hypId]++;
        hypId++;
      }
    }
    std::cout << std::endl;

    // save the result to std::map<nodeId, blockid>
    std::map<mt_kahypar_hypernode_id_t, mt_kahypar_partition_id_t> node_block_io;
    std::map<mt_kahypar_hypernode_id_t, mt_kahypar_partition_id_t> node_block;
    for ( auto& item : vBoundaries )
    {
      mt_kahypar_hypernode_id_t hypId = 0;
      while ( hypId < item.second.size() )
      {
        if ( node_block_io.count( item.second[hypId] ) > 0 )
        {
          hypId++;
          continue;
        }
        else
        {
          node_block_io[item.second[hypId] - 1] = mt_kahypar_block_id( partitioned_hg, item.second[hypId] - 1 );
          hypId++;
        }
      }
    }
    std::cout << "Size of all the boundary nodes is " << node_block_io.size() << std::endl;
    node_block = convertToMap( partition, mt_kahypar_num_hypernodes( hypergraph ) );

    std::cout << "Original aig PI num: " << aig.num_pis() << std::endl;
    std::cout << "Original aig PO num: " << aig.num_pos() << std::endl;

    auto vAigs = aig_p.construct_from_partition( 2, vBoundaries, node_block_io, node_block );

    // Now insert all the aigs back to original ntk and check the equivalence
    for ( auto& aig_part : vAigs )
    {

      klut_network klut_a = lut_map( std::get<0>( aig_part ) );
      mig_network mig_IR;
      convert_klut_to_graph<mig_network>( mig_IR, klut_a );
      klut_network klut_m = lut_map( mig_IR );
      aig_network aig_p_t_b;
      convert_klut_to_graph<aig_network>( aig_p_t_b, klut_m );
      std::get<0>( aig_part ) = aig_p_t_b;

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

    auto aig_clear = cleanup_dangling( aig );
    color_view f_c_aig{ aig_clear };
    assert( count_reachable_dead_nodes( f_c_aig ) == 0u );
    assert( network_is_acyclic( f_c_aig ) );
    aig = aig_clear;

    auto final_gate_num = aig.num_gates();
    std::cout << "Original gate number " << ori_gate_num << " Final gate number " << final_gate_num << std::endl;

    // equivalence check
    const auto cec1 = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    assert( cec1 == true );
    std::cout << "*****EQ CHEKCED*****" << std::endl;

    // Extract Block Weights
    auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>( 2 );
    mt_kahypar_get_block_weights( partitioned_hg, block_weights.get() );

    // Compute Metrics
    const double imbalance = mt_kahypar_imbalance( partitioned_hg, context );
    const int km1 = mt_kahypar_km1( partitioned_hg );

    // Output Results
    std::cout << "====================\nPartitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;
    // std::cout << "Weight of Block 2 = " << block_weights[2] << std::endl;
    std::cout << std::endl;

    // You cannot give a directory that doesn't exist here, std::ofstream doesn't create new directory by default
    mt_kahypar_status_t status_out = mt_kahypar_write_partition_to_file( partitioned_hg,
                                                                         fmt::format( "{}/{}hypOut", pHyOutS, benchmark ).c_str(),
                                                                         &error );
    assert( status_out == SUCCESS );

    // alternative without file write out, using:
    // mt_kahypar_partitioned_hypergraph_t partitioned_hg_sep = mt_kahypar_create_partitioned_hypergraph(hypergraph, context, 3, partition.get(), &error);

    // if ( remove( fmt::format( "{}/tmp.hmetis", pHyOutS ).c_str() ) == 0 )
    //{
    //   std::cout << "Successfully deleted the tmp file." << std::endl;
    // }
    // else
    //{
    //   std::cout << "Failed to delete the tmp file." << std::endl;
    // }

    mt_kahypar_free_context( context );
    mt_kahypar_free_hypergraph( hypergraph );
    mt_kahypar_free_partitioned_hypergraph( partitioned_hg );
  }
  return 0;
}