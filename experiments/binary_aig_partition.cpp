/*!
  \file binary_aig_partition.cpp
  \brief Read a binary AIG file, partition it using mt-KaHyPar, and write partition AIGs to disk.

  Usage: ./binary_aig_partition <input.aig> <num_blocks> <output_dir>

  Each partition is saved as:
    - <output_dir>/partition_N.aig        (binary AIG)
    - <output_dir>/partition_N_meta.json  (metadata for stitching)

  \author Jingren Wang
*/

#include "experiments.hpp"
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <lorina/aiger.hpp>
#include <memory>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/partition_view.hpp>
#include <mtkahypar.h>
#include <nlohmann/json.hpp>
#include <string>
#include <thread>
#include <vector>

using namespace mockturtle;

int main( int argc, char* argv[] )
{
  if ( argc != 4 )
  {
    fmt::print( "Usage: {} <input.aig> <num_blocks> <output_dir>\n", argv[0] );
    return 1;
  }

  std::string input_file = argv[1];
  int num_blocks = std::stoi( argv[2] );
  std::string output_dir = argv[3];

  std::filesystem::create_directories( output_dir );

  aig_network aig;
  if ( lorina::read_aiger( input_file, aiger_reader( aig ) ) != lorina::return_code::success )
  {
    fmt::print( "[e] Failed to read AIG file: {}\n", input_file );
    return 1;
  }

  fmt::print( "[i] Read AIG: {} PIs, {} POs, {} gates\n",
              aig.num_pis(), aig.num_pos(), aig.num_gates() );

  std::string hmetis_file = fmt::format( "{}/_tmp_partition.hmetis", output_dir );

  partition_view_params ps;
  ps.file_name = hmetis_file;
  ps.num_blocks = num_blocks;

  partition_view aig_p{ aig, ps };

  mt_kahypar_error_t error{};
  mt_kahypar_initialize(
      std::thread::hardware_concurrency(),
      true );

  mt_kahypar_context_t* context = mt_kahypar_context_from_preset( DETERMINISTIC );
  mt_kahypar_set_partitioning_parameters( context,
                                          num_blocks, 0.03,
                                          KM1 );
  mt_kahypar_set_seed( 42 );
  mt_kahypar_set_context_parameter( context, VERBOSE, "0", &error );

  mt_kahypar_hypergraph_t hypergraph =
      mt_kahypar_read_hypergraph_from_file( hmetis_file.c_str(),
                                            context, HMETIS, &error );
  if ( hypergraph.hypergraph == nullptr )
  {
    fmt::print( "[e] Failed to load hypergraph: {}\n", error.msg );
    return 1;
  }

  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_partition( hypergraph, context, &error );
  if ( partitioned_hg.partitioned_hg == nullptr )
  {
    fmt::print( "[e] Partitioning failed: {}\n", error.msg );
    return 1;
  }

  auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(
      mt_kahypar_num_hypernodes( hypergraph ) );
  mt_kahypar_get_partition( partitioned_hg, partition.get() );

  auto vAigs = aig_p.construct_from_partition( num_blocks, partition, hypergraph );

  for ( int i = 0; i < num_blocks; i++ )
  {
    auto& aig_part = vAigs[i];
    auto& sub_aig = std::get<0>( aig_part );
    auto& inputs = std::get<1>( aig_part );
    auto& outputs = std::get<2>( aig_part );
    auto& gates = std::get<3>( aig_part );

    std::string part_file = fmt::format( "{}/partition_{}.aig", output_dir, i );
    write_aiger( sub_aig, part_file );

    nlohmann::json meta;
    meta["partition_index"] = i;
    meta["partition_file"] = fmt::format( "partition_{}.aig", i );
    meta["num_gates"] = sub_aig.num_gates();
    meta["num_pis"] = sub_aig.num_pis();
    meta["num_pos"] = sub_aig.num_pos();

    for ( auto const& n : inputs )
    {
      meta["input_node_indices"].push_back(
          static_cast<uint64_t>( aig.node_to_index( n ) ) );
    }
    for ( auto const& s : outputs )
    {
      meta["output_signal_indices"].push_back(
          static_cast<uint64_t>( aig.node_to_index( aig.get_node( s ) ) ) );
      meta["output_signal_complemented"].push_back(
          aig.is_complemented( s ) );
    }
    for ( auto const& n : gates )
    {
      meta["gate_node_indices"].push_back(
          static_cast<uint64_t>( aig.node_to_index( n ) ) );
    }

    std::string meta_file = fmt::format( "{}/partition_{}_meta.json", output_dir, i );
    std::ofstream ofs( meta_file );
    ofs << meta.dump( 2 ) << "\n";
    ofs.close();

    fmt::print( "[i] Partition {}: {} gates, {} PIs, {} POs -> {}\n",
                i, sub_aig.num_gates(), sub_aig.num_pis(), sub_aig.num_pos(),
                part_file );
  }

  std::filesystem::remove( hmetis_file );

  mt_kahypar_free_context( context );
  mt_kahypar_free_hypergraph( hypergraph );
  mt_kahypar_free_partitioned_hypergraph( partitioned_hg );

  fmt::print( "[i] Partitioning complete. {} partitions saved to {}\n",
              num_blocks, output_dir );
  return 0;
}
