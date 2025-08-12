# SRonGIG: Split-and-Reunion on GIG

<p align="center">
   <picture>
     <img src="./logo/SRONGIG.jpeg" width="50%" alt="SRonGIG logo">
   </picture>
</p>

Split-and-Reunion on Graph-Inverter-Graph is an interface we implemented to enable partition process in logic synthesis.

It is based on [mockturtle](https://github.com/lsils/mockturtle) and will finally be merged into a self-maintained version.

## What & Why & how
It is a graph based interface that could easily be used with other mockturtle interfaces. We implement it to assist future work which based on graph partitioning. We use SOTA [mt-KaHypa](https://github.com/kahypar/mt-kahypar) as bridge. The efficient data format representing the hypergraph is called [hMetis format](https://course.ece.cmu.edu/~ee760/760docs/hMetisManual.pdf) (it is highly recommended to directly read Figure 5 and related paragraph to get an intuitive idea).

## Usage
### Dump the hMetis file format
Now you could simply include the partition view to dump the hMetis file format. Eg,
```cpp
#include <mockturtle/views/partition_view.hpp>
aig_network aig;
...
partition_view_params ps;
ps.file_name = fmt::format( "{}/test.hmetis", PARTITION_TEST_PATH );
partition_view aig_p{ aig, ps };
```
It will dump a hMetis file without weight on edge or vertices by default, you could enable the simple model by turn it/them on. Eg,
```cpp
#include <mockturtle/views/partition_view.hpp>
aig_network aig;
...
ps.si_w_on_hyperedges = true;
ps.si_w_on_vertices = true;
ps.file_name = fmt::format( "{}/test_edge_vertices_weight.hmetis", PARTITION_TEST_PATH );
partition_view aig_p_e_v_w{ aig, ps };
```
Note that the weight should be modified by the target of the problem, this is just a showcase that the hyperedges have the weight of the fanout number of each source node, and vertices have the weight of 1.

> [!NOTE]
> We support simple weights on vertices and edges, but users should customize them for their own need.

### Partition process
> [!TIP]
> Since mt-KaHypa is used in this view, read the experiment of [experiments/read_partition_data_format.cpp](experiments/read_partition_data_format.cpp) is strongly recommended since it also includes interface examples of mt-KaHypa, and give you a minimum example of how to do the partition.

#### Get the partitions
Use the partition is quite handy now, after you provide the number of blocks, `partition` and `hypergraph`, the AIG partition information will be saved in a vector,
```cpp
auto vAigs = aig_p.construct_from_partition( ps.num_blocks, partition, hypergraph );
```
which has the type of
```cpp
std::vector<std::tuple<aig_network, std::vector<node>, std::vector<signal>, std::vector<node>>>
```
and the pure AIG can be fetched from each one of the tuple by `std::get<0>( aig_part )`.

#### Stitch back to original AIG
After dealing with each part of the partition, you could insert it back with a simple
```cpp
aig_p.insert_back( aig_part );
```
It is **vital** that you remove the **dangling nodes** after all partitions insert back to original AIG, since it will influence the equivalence results.

## Build Experiment
Please follow the instruction in mockturtle, we are not breaking any mockturtle project structure. We give a showcase in [experiments/reader_simple_partition.cpp](experiments/reader_simple_partition.cpp). So a simple build and test on partition would be:

```bash
mkdir build
cd build
cmake -DMOCKTURTLE_EXPERIMENTS=ON -DCMAKE_BUILD_TYPE=Debug -DMOCKTURTLE_TEST=ON ..
make reader_simple_partition
./experiments/reader_simple_partition
```

## Automatic Latex Data Collection
The experimental framework leverages Mockturtle's robust JSON generation capabilities for data output. We have implemented a Python-based processing pipeline that automatically transforms this collected data into LaTeX-compatible formats (see jupyternotebook [experiments/data_collect/DataCollectionToLatex.ipynb](experiments/data_collect/DataCollectionToLatex.ipynb) for implementation details).

> [!NOTE]
> This is a simple version of data transformation, for detailed interfaces, check [pandas.DataFrame.to_latex](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_latex.html).