/*!
  \file binary_aig_stitch.cpp
  \brief Read previously partitioned binary AIG files and stitch them back into a whole AIG.

  Usage: ./binary_aig_stitch <original.aig> <partition_dir> <output.aig>

  Expects the partition directory to contain:
    - partition_N.aig        (binary AIG files, from binary_aig_partition)
    - partition_N_meta.json  (metadata files, from binary_aig_partition)

  It is vital to remove dangling nodes after stitching. This tool does it
  automatically (as noted in the README).

  \author Jingren Wang
*/

#include "experiments.hpp"
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/utils/network_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using namespace mockturtle;

int main( int argc, char* argv[] )
{
  if ( argc != 4 )
  {
    fmt::print( "Usage: {} <original.aig> <partition_dir> <output.aig>\n", argv[0] );
    return 1;
  }

  std::string original_file = argv[1];
  std::string partition_dir = argv[2];
  std::string output_file = argv[3];

  aig_network aig;
  if ( lorina::read_aiger( original_file, aiger_reader( aig ) ) != lorina::return_code::success )
  {
    fmt::print( "[e] Failed to read original AIG: {}\n", original_file );
    return 1;
  }

  fmt::print( "[i] Original AIG: {} PIs, {} POs, {} gates\n",
              aig.num_pis(), aig.num_pos(), aig.num_gates() );

  std::vector<std::string> meta_files;
  for ( int i = 0;; i++ )
  {
    std::string meta_file = fmt::format( "{}/partition_{}_meta.json", partition_dir, i );
    if ( !std::filesystem::exists( meta_file ) )
    {
      break;
    }
    meta_files.push_back( meta_file );
  }

  if ( meta_files.empty() )
  {
    fmt::print( "[e] No partition metadata files found in {}\n", partition_dir );
    return 1;
  }

  fmt::print( "[i] Found {} partitions\n", meta_files.size() );

  for ( auto const& meta_file : meta_files )
  {
    std::ifstream ifs( meta_file );
    nlohmann::json meta = nlohmann::json::parse( ifs );
    ifs.close();

    std::string part_file = fmt::format( "{}/{}", partition_dir,
                                         meta["partition_file"].get<std::string>() );

    aig_network part_aig;
    if ( lorina::read_aiger( part_file, aiger_reader( part_aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[e] Failed to read partition AIG: {}\n", part_file );
      return 1;
    }

    std::vector<aig_network::signal> i_sigs;
    for ( auto const& idx : meta["input_node_indices"] )
    {
      auto node_idx = static_cast<uint32_t>( idx.get<uint64_t>() );
      auto n = aig.index_to_node( node_idx );
      i_sigs.push_back( aig.make_signal( n ) );
    }

    std::vector<aig_network::signal> old_outputs;
    auto const& sig_indices = meta["output_signal_indices"];
    auto const& sig_compl = meta["output_signal_complemented"];
    for ( size_t j = 0; j < sig_indices.size(); j++ )
    {
      auto node_idx = static_cast<uint32_t>( sig_indices[j].get<uint64_t>() );
      bool complemented = sig_compl[j].get<bool>();
      auto n = aig.index_to_node( node_idx );
      auto s = complemented ? !aig.make_signal( n ) : aig.make_signal( n );
      old_outputs.push_back( s );
    }

    color_view c_aig{ aig };
    assert( count_reachable_dead_nodes( c_aig ) == 0u );

    uint32_t counter = 0u;
    insert_ntk( aig, i_sigs.begin(), i_sigs.end(), part_aig,
                [&]( aig_network::signal const& _new ) {
                  assert( !c_aig.is_dead( c_aig.get_node( _new ) ) );
                  auto const _old = old_outputs.at( counter++ );
                  if ( _old == _new )
                  {
                    return true;
                  }

                  if ( _old != _new )
                  {
                    aig.substitute_node( aig.get_node( _old ),
                                         aig.is_complemented( _old ) ? !_new : _new );
                  }
                  return true;
                } );

    fmt::print( "[i] Stitched partition {} ({} gates)\n",
                meta["partition_index"].get<int>(), part_aig.num_gates() );
  }

  auto aig_clear = cleanup_dangling( aig );
  color_view f_c_aig{ aig_clear };
  assert( count_reachable_dead_nodes( f_c_aig ) == 0u );
  assert( network_is_acyclic( f_c_aig ) );
  aig = aig_clear;

  write_aiger( aig, output_file );

  fmt::print( "[i] Stitched AIG: {} PIs, {} POs, {} gates -> {}\n",
              aig.num_pis(), aig.num_pos(), aig.num_gates(), output_file );
  return 0;
}
