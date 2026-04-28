/*!
  \file binary_aig_partition.cpp
  \brief Read a binary AIG file, partition it using mt-KaHyPar, and write partition AIGs to disk.

  Usage: ./binary_aig_partition <input.aig> <num_blocks> <output_dir> [options]

  Options:
    --seed <int>           Random seed (default: 42)
    --epsilon <float>      Imbalance tolerance (default: 0.03)
    --objective <str>      Objective: cut, km1, soed (default: km1)
    --preset <str>         Preset: deterministic, default, quality, highest_quality, large_k (default: deterministic)
    --vcycles <int>        Number of V-cycles (default: 0)
    --edge-weight          Enable hyperedge weights (fanout-based)
    --vertex-weight         Enable vertex weights (uniform 1)

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

mt_kahypar_objective_t parse_objective( std::string const& s )
{
  if ( s == "cut" )
    return CUT;
  if ( s == "soed" )
    return SOED;
  return KM1;
}

mt_kahypar_preset_type_t parse_preset( std::string const& s )
{
  if ( s == "default" )
    return DEFAULT;
  if ( s == "quality" )
    return QUALITY;
  if ( s == "highest_quality" )
    return HIGHEST_QUALITY;
  if ( s == "large_k" )
    return LARGE_K;
  return DETERMINISTIC;
}

int main( int argc, char* argv[] )
{
  if ( argc < 4 )
  {
    fmt::print( "Usage: {} <input.aig> <num_blocks> <output_dir> [options]\n", argv[0] );
    fmt::print( "Options:\n" );
    fmt::print( "  --seed <int>           Random seed (default: 42)\n" );
    fmt::print( "  --epsilon <float>      Imbalance tolerance (default: 0.03)\n" );
    fmt::print( "  --objective <str>      cut, km1, soed (default: km1)\n" );
    fmt::print( "  --preset <str>         deterministic, default, quality, highest_quality, large_k (default: deterministic)\n" );
    fmt::print( "  --vcycles <int>        Number of V-cycles (default: 0)\n" );
    fmt::print( "  --edge-weight          Enable hyperedge weights\n" );
    fmt::print( "  --vertex-weight        Enable vertex weights\n" );
    return 1;
  }

  std::string input_file = argv[1];
  int num_blocks = std::stoi( argv[2] );
  std::string output_dir = argv[3];

  size_t seed = 42;
  double epsilon = 0.03;
  mt_kahypar_objective_t objective = KM1;
  mt_kahypar_preset_type_t preset = DETERMINISTIC;
  int vcycles = 0;
  bool edge_weight = false;
  bool vertex_weight = false;

  for ( int i = 4; i < argc; i++ )
  {
    std::string arg = argv[i];
    if ( arg == "--seed" && i + 1 < argc )
    {
      seed = std::stoull( argv[++i] );
    }
    else if ( arg == "--epsilon" && i + 1 < argc )
    {
      epsilon = std::stod( argv[++i] );
    }
    else if ( arg == "--objective" && i + 1 < argc )
    {
      objective = parse_objective( argv[++i] );
    }
    else if ( arg == "--preset" && i + 1 < argc )
    {
      preset = parse_preset( argv[++i] );
    }
    else if ( arg == "--vcycles" && i + 1 < argc )
    {
      vcycles = std::stoi( argv[++i] );
    }
    else if ( arg == "--edge-weight" )
    {
      edge_weight = true;
    }
    else if ( arg == "--vertex-weight" )
    {
      vertex_weight = true;
    }
    else
    {
      fmt::print( "[w] Unknown option: {}\n", arg );
    }
  }

  fmt::print( "[i] Parameters: seed={}, epsilon={}, objective={}, preset={}, vcycles={}, edge_weight={}, vertex_weight={}\n",
              seed, epsilon,
              objective == CUT ? "cut" : ( objective == SOED ? "soed" : "km1" ),
              preset == DETERMINISTIC ? "deterministic" : ( preset == DEFAULT ? "default" : ( preset == QUALITY ? "quality" : ( preset == HIGHEST_QUALITY ? "highest_quality" : "large_k" ) ) ),
              vcycles, edge_weight, vertex_weight );

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
  ps.seed = seed;
  ps.epsilon = epsilon;
  ps.si_w_on_hyperedges = edge_weight;
  ps.si_w_on_vertices = vertex_weight;

  partition_view aig_p{ aig, ps };

  mt_kahypar_error_t error{};
  mt_kahypar_initialize(
      std::thread::hardware_concurrency(),
      true );

  mt_kahypar_context_t* context = mt_kahypar_context_from_preset( preset );
  mt_kahypar_set_partitioning_parameters( context,
                                          num_blocks, epsilon,
                                          objective );
  mt_kahypar_set_seed( seed );
  mt_kahypar_set_context_parameter( context, VERBOSE, "0", &error );

  if ( vcycles > 0 )
  {
    mt_kahypar_set_context_parameter( context, NUM_VCYCLES,
                                      std::to_string( vcycles ).c_str(), &error );
  }

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

  nlohmann::json partition_info;
  partition_info["input_file"] = input_file;
  partition_info["num_blocks"] = num_blocks;
  partition_info["seed"] = seed;
  partition_info["epsilon"] = epsilon;
  partition_info["objective"] = objective == CUT ? "cut" : ( objective == SOED ? "soed" : "km1" );
  partition_info["preset"] = preset == DETERMINISTIC ? "deterministic" : ( preset == DEFAULT ? "default" : ( preset == QUALITY ? "quality" : ( preset == HIGHEST_QUALITY ? "highest_quality" : "large_k" ) ) );
  partition_info["vcycles"] = vcycles;
  partition_info["edge_weight"] = edge_weight;
  partition_info["vertex_weight"] = vertex_weight;

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

    partition_info["partitions"].push_back( { { "index", i },
                                              { "gates", sub_aig.num_gates() },
                                              { "pis", sub_aig.num_pis() },
                                              { "pos", sub_aig.num_pos() } } );
  }

  {
    std::string info_file = fmt::format( "{}/partition_info.json", output_dir );
    std::ofstream ofs( info_file );
    ofs << partition_info.dump( 2 ) << "\n";
    ofs.close();
  }

  std::filesystem::remove( hmetis_file );

  auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>( num_blocks );
  mt_kahypar_get_block_weights( partitioned_hg, block_weights.get() );
  const double imbalance = mt_kahypar_imbalance( partitioned_hg, context );
  const int km1 = mt_kahypar_km1( partitioned_hg );
  const int cut = mt_kahypar_cut( partitioned_hg );
  const int soed = mt_kahypar_soed( partitioned_hg );

  fmt::print( "[i] Partitioning metrics: imbalance={:.4f}, km1={}, cut={}, soed={}\n",
              imbalance, km1, cut, soed );
  for ( int i = 0; i < num_blocks; i++ )
  {
    fmt::print( "[i]   Block {}: weight={}\n", i, block_weights[i] );
  }

  mt_kahypar_free_context( context );
  mt_kahypar_free_hypergraph( hypergraph );
  mt_kahypar_free_partitioned_hypergraph( partitioned_hg );

  fmt::print( "[i] Partitioning complete. {} partitions saved to {}\n",
              num_blocks, output_dir );
  return 0;
}
