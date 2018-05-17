# Eigenvector-based identification of bipartite subgraphs
## Summary
* The script BIPARTITE_SUBGRAPH_IDENTIFICATION returns the relative size of bipartition (the ratio between the number of edges in the bipartition to the total number of edges in the orginal graph) for a simple, connected, and undirected graph using: 
   
   * A combinatorial method following Erdös bound
   
   * Sign-based partioning of eigenvectors belongs to the 
   
      *  smallest eigenvalues of the adjacency and signless Laplacian matrix
      *  largest eigenvalues of the Laplacian and normalized Laplacian matrix
      
   * Estrada & Gómez-Gardeñes  edge bipartivity [Estrada, E. and Gómez-Gardeñes, J., Physica D: Nonlinear Phenomena, 323, pp.57-63, 2016]
   
   * Two new measures based on the adjancency and the normalized laplacian matrix [Paul, D. and Stevanovic, D., (2018) (in preparation)] 

## What is a bipartite graph?
