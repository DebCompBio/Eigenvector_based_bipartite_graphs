# Eigenvector-based identification of bipartite subgraphs


## What is a bipartite graph?
A graph *G*(*V*, *E*) with a set of vertices *V* and a set of edges *E*
is bipartite if there exists a partition *V* = *X* ∪ *Y*, *X* ∩ *Y* = ∅
such that every edge *e* ∈ *E* has one end in *X* and another in *Y*. It
is a classical result that a graph is bipartite if and only if it does
not contain a cycle of odd length as a subgraph.


## Description of the script
* The script BIPARTITE_SUBGRAPH_IDENTIFICATION returns the relative size of bipartition (the ratio between the number of edges in the bipartition to the total number of edges in the orginal graph) for a simple, connected, and undirected graph using: 
   
   * A combinatorial method following Erdös bound (0.5)[1] 
   
   * Sign-based partitioning of eigenvectors belong to the 
   
      *  smallest eigenvalues of the adjacency and signless Laplacian matrix
      *  largest eigenvalues of the Laplacian and normalized Laplacian matrix
      
   *  Edge bipartivity by Estrada & Gómez-Gardeñes [2]
   
   * Two new measures based on the adjancency and the normalized laplacian matrix [3] 

[1] Erdös P. On some extremal problems in graph theory. Israel Journal of Mathematics 1965;3(2):113–6.

[2] Estrada, E. and Gómez-Gardeñes, J., Physica D: Nonlinear Phenomena, 323, pp.57-63, 2016

[3] Paul, D. and Stevanovic, D., (2018) (in preparation)

## Dependencies

* networkx (https://networkx.github.io/documentation/networkx-1.10/index.html)
* numpy
* random
