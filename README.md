# Eigenvector-based identification of bipartite subgraphs


## What is a bipartite graph?
A graph *G*(*V*, *E*) with a set of vertices *V* and a set of edges *E*
is bipartite if there exists a partition *V* = *X* ∪ *Y*, *X* ∩ *Y* = ∅
such that every edge *e* ∈ *E* has one end in *X* and another in *Y*. It
is a classical result that a graph is bipartite if and only if it does
not contain a cycle of odd length as a subgraph.


## Description of the script
* The script BIPARTITE_SUBGRAPH_IDENTIFICATION returns the relative size of bipartition (the ratio between the number of edges in the bipartition to the total number of edges in the orginal graph) for a simple, connected, and undirected graph using: 
   
   * A combinatorial method (local switching search algorithm - [1]) following Erdös bound (0.5)[2] 
   
   * Sign-based partitioning of eigenvectors belong to the 
   
      *  smallest eigenvalues of the adjacency and signless Laplacian matrix
      *  largest eigenvalues of the Laplacian and normalized Laplacian matrix
      
   *  Edge bipartivity by Estrada & Gómez-Gardeñes [3]
   
   * Two new measures based on the adjancency and the normalized laplacian matrix [4] 

[1] Bylka, S., Idzik, A. and Tuza, Z.,Discrete Mathematics,  1999;194(1-3), pp.39-58.

[2] Erdös P., Israel Journal of Mathematics 1965;3(2):113–6.

[3] Estrada, E. and Gómez-Gardeñes, J., Physica D: Nonlinear Phenomena, 2016;323, pp.57-63, 

[4] Paul, D., & Stevanović, D. (2019). Eigenvector-based identification of bipartite subgraphs. Discrete Applied Mathematics.
 

## How to run?
* Input: 
  ```shell 
  python bipartite_subgraph_identification.py -h 
  ```
* Output: 
```shell
usage: bipartite_subgraph_identification.py [-h] [--g [G]] [--p_er P_ER P_ER]
                                            [--p_ba P_BA P_BA]
                                            [--p_rg P_RG P_RG]
                                            [--p_ws P_WS P_WS P_WS] [--n N]
                                            [--ng NG] [--o O]

optional arguments:
  -h, --help            show this help message and exit
  --g [G], --graph [G]  graph models: ER, BA, RG, WS
  --p_er P_ER P_ER, --param_er P_ER P_ER
                        parameter for ER graph model: min max
  --p_ba P_BA P_BA, --param_ba P_BA P_BA
                        parameter for BA graph model: min max
  --p_rg P_RG P_RG, --param_rg P_RG P_RG
                        parameter for RG graph model: min max
  --p_ws P_WS P_WS P_WS, --param_ws P_WS P_WS P_WS
                        parameters k (nearest neighbors) & r (probability of
                        rewiring)for WS graph model
  --n N, --nodes N      number of nodes
  --ng NG, --nog NG     number of graphs
  --o O, --output O     name of the output folder
```
* Example command: 
```shell 
   python bipartite_subgraph_identification.py --g ER --p_er 0.2 1.0 --n 20 --ng 2 --o test_output 
 ```
The above command will create the following 8 files inside the folder <./test_output/graph_type/> for 2 ER type graphs each with 20 vertices:

  * ratioErdos.txt    : ratio obtained from the combinatorial method
  * ratioA.txt        : ratio obtained from the adjacency eigenvectors
  * ratioQ.txt        : ratio obtained from the signless Laplacian eigenvectors
  * ratioL.txt        : ratio obtained from the Laplacian eigenvectors
  * ratioNL.txt       : ratio obtained from the normalized Laplacian eigenvectors 
  * ratioEstrada.txt  : ratio obtained from the Estrada-Gómez-Gardeñes edge bipartivity measure
  * ratioAnew.txt     : ratio obtained from the proposed edge bipartivity measure based on the adjacency matrix
  * ratioNLnew.txt    : ratio obtained from the proposed edge bipartivity measure based on the normalized laplacian matrix

## Dependencies

* networkx (https://networkx.github.io/documentation/networkx-1.10/index.html)
* numpy
* random
* matplotlib (for plotting)
