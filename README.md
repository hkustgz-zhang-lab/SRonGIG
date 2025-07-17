# SRonGIG: Split-and-Reunion on GIG

Split-and-Reunion on Graph-Inverter-Graph is an interface we implemented to enable partition process in logic synthesis.

It is based on [mockturtle](https://github.com/lsils/mockturtle) and will finally be merged into a self-maintained version.

## What & Why & how
It is a graph based interface that could easily be used with other mockturtle interfaces. We implement it to assist future work which based on graph partitioning. We use SOTA [mt-KaHypa](https://github.com/kahypar/mt-kahypar) as bridge. The efficient data format representing the hypergraph is called [hMetis format](https://course.ece.cmu.edu/~ee760/760docs/hMetisManual.pdf) (it is highly recommended to directly read Figure 5 and related paragraph to get an intuitive idea).

## Build
Please follow the instruction in mockturtle, we are not breaking any mockturtle project structure. We give a showcase in [experiments/reader_simple_partition.cpp](experiments/reader_simple_partition.cpp). So a simple build and test on partition would be:

```bash
mkdir build
cd build
cmake -DMOCKTURTLE_EXPERIMENTS=ON -DCMAKE_BUILD_TYPE=Debug -DMOCKTURTLE_TEST=ON ..
make reader_simple_partition
./experiments/reader_simple_partition
```
