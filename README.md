# Maximum-kDC: Maximum k-Defective Clique Computation

This repository implements the maximum k-defective clique computation algorithm proposed in our SIGMOD 2024 paper. If you are using the code, please cite our paper.
<pre>
Lijun Chang.
<a href="https://lijunchang.github.io/pdf/2024-Maximum-kDC.pdf">Efficient Maximum k-Defective Clique Computation with Improved Time Complexity.</a>
SIGMOD (2024)
</pre>

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "kDC", which corresponds to the kDC algorithm.

## Run the code

```sh
$ ./kDC -g {path_to_graph} -k {k_value}
```

An example of computing the exact maximum 3-defective clique for the dataset CA-GrQc is as follows
```sh
$ ./kDC -g datasets/CA-GrQc -k 3
```

## Data format
Two data formats are supported. The default data format is "edges.txt", which contains a list of undirected edges represented as vertex pairs. The first line contains two numbers n and m, representing the number of vertices and the number of undirected edges, respectively. Note that, the vertex ids must be between 0 and n-1.

The more time-efficient format is the binary format; to read the input graph from this format, please add "-b" when running the code. Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://lijunchang.github.io/Cohesive_subgraph_book/datasets](https://lijunchang.github.io/Cohesive_subgraph_book/datasets)

## Get datasets
Real-world graphs collection: [http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz]

Facebook graphs collection: [https://networkrepository.com/socfb.php]
