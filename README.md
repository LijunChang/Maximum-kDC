# Maximum_kDC: Maximum k-Defective Clique Computation

This repository implements the maximum k-defective clique computation algorithm proposed in our SIGMOD 2024 paper. 

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

## Get other datasets
Real-world graphs collection: http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz
Facebook graphs collection: https://networkrepository.com/socfb.php
