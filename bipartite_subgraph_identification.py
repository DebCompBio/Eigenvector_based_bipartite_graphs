#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 15:35:41 2018

@author: Paul, Debdas 

"""
""" Description 
    
    The script COMBINATORIAL_EIGENVECTOR_APPROACH contains a list of routines 
    to obtail large bipartite subgraphs of a graph based on 
    
    1. Combinatorial approach: 
        
       [Ref:] Erdös P. On some extremal problems in graph theory, 
                 Israel Journal of Matematics,  (1965); 3(2):pp. 113-6
    2. Sign-pattern of eigenvector based approach 
    
    3. Estarda's edge bipartivity:
        
       [Ref:] Estrada, E., & Gómez-Gardeñes, J. (2016). Network bipartivity and the transportation efficiency of European passenger airlines. 
              Physica D: Nonlinear Phenomena, 323, 57-63.
    4. Eigenevector based two new approaches mentioned in:
        
       [Ref:] Paul, D. & Stevanovic, D. (2018). Eigenvector-based identification of 
               bipartite subgraphs, (submitted)
    
    
    INPUT :  A simple, connected, and undirected graph G
    
    OUTPUT: A bipartiiton of vertices with the ratio between the edges in the bipartition
            to the total number of edges in the graph G
    
"""


    

""" 
    Graph models and paramters

    1. Erdös-Renyi
    2. Barabasi-Albert
    3. Random-Geometric
    4. Watts-Strogatz
    
    
    Parameter specifications:
    
    1.       
        n (int) – The number of nodes.
        p (float) – Probability for edge creation. should be > ln(n)/n to be connected
        seed (int, optional) – Seed for random number generator (default=None).
        directed (bool, optional (default=False)) – If True, this function returns a directed 

    2. 	
        n (int) – Number of nodes
        m (int) – Number of edges to attach from a new node to existing nodes
        seed (int, optional) – Seed for random number generator (default=None). 
    
    3.  
    
        n (int or iterable) – Number of nodes or iterable of nodes
        radius (float) – Distance threshold value
        dim (int, optional) – Dimension of graph
        pos (dict, optional) – A dictionary keyed by node with node positions as values.
        p (float, optional) – Which Minkowski distance metric to use. p has to meet the condition 1 <= p <= infinity. 
   
    
    4. 
       n (int) – The number of nodes
       k (int) – Each node is joined with its k nearest neighbors in a ring topology.
       p (float) – The probability of rewiring each edge
       tries (int) – Number of attempts to generate a connected graph.
       seed (int, optional) – The seed for random number generator.
"""





""" import packages """

import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import argparse
import os

def make_partition(lst,n):
    """ Make random partition of vertices """
    
    division = len(lst) / float(n) 
    return [lst[int(round(division * i)):
            int(round(division * (i + 1)))] for i in range(n)] 

def make_shuffle(lst):  

    return random.sample(lst,len(lst));      
    
def update_partition(vertex_u,node_partition,choice_of_partition):

    """ Move a vertex from one partition to another"""
    node_partition[choice_of_partition].remove(vertex_u)
    node_partition[1-int(choice_of_partition)].extend([vertex_u])
    return node_partition
    
    
def global_vertex_list(G):
   track_movement_list={};
   
   for node in  list_of_nodes(G):
       
       track_movement_list[node]=0;
       
   return track_movement_list 

def update_global_list(vertex_u,track_movement_list):
    
    """ Update the global list of vertices to check
        if there will be further movement of vertices"""
        
    track_movement_list[vertex_u]=1;
    return track_movement_list;
       
def  adjacency_matrix(G):
    
    return nx.to_numpy_matrix(G);

def  laplacian_matrix(G):
    
      return np.matrix(degree_matrix(G) - adjacency_matrix(G)) 

def signless_laplacian_matrix(G):

      return np.matrix(degree_matrix(G) + adjacency_matrix(G))

def normalized_laplacian_matrix(G):
    
     diag_sqrt = degree_matrix(G);
     
     for row in range(len(diag_sqrt)):
         
         diag_sqrt[row,row] = pow(diag_sqrt[row,row],-0.5);
         
     return np.matrix(diag_sqrt*laplacian_matrix(G)*diag_sqrt)    
 
def degree_matrix(G):
     transposed_adjacency_matrix=compute_transpose(adjacency_matrix(G))
     identity_matrix=return_identity_matrix(len(G.nodes()))
     
     for diag_elements in range(len(G.nodes())):
         
         identity_matrix[diag_elements,diag_elements]=transposed_adjacency_matrix[0,diag_elements];
     
     return identity_matrix
         

def list_of_nodes(G):
    
    """ Return the list of nodes of graph G """
    
    return list(G.nodes)

def compute_transpose(matrix):
    
    return np.transpose(np.sum(matrix,1))

def return_identity_matrix(length_of_list):
    
    return np.eye(length_of_list)


def compute_edges_within_partition(partition,adj_matrix):
    """ Compute the neighbors of a vertex u within its member 
        partion"""
    
    extracted_adj_matrix=adj_matrix[np.ix_(partition, partition)] 
    
    return np.sum(extracted_adj_matrix,axis=1)
    
    
def compute_edges_across_partition(partition,adj_matrix):
    """ Compute the neighbors of a vertex u across its member 
        partion"""
    return compute_edges(adj_matrix)[np.ix_(partition, [0])]-compute_edges_within_partition(partition,adj_matrix);

def compute_edges(Adj):
    
    return np.sum(Adj,axis=1)


def check_connected(G):
    
    """Check if the graph is connected or not """
    
    return nx.is_connected(G);
    
def edwards_ratio(G):
    
   return float(G.number_of_edges())/float(2) + (float(len(G.nodes()))-1.0)/4.0;
    

""" Graph generation routines """

def generate_er(n,p):
    
    gr=nx.erdos_renyi_graph(n, p, seed=None, directed=False);
    print ("Graph type: Erdös-Renyi")
    print ("\n number of nodes:", nx.number_of_nodes(gr))
    print ("\n number of edges:", nx.number_of_edges(gr))
    return gr

def generate_ba(n,m):
      
    return  nx.barabasi_albert_graph(n, m, seed=None);
    
def generate_rg(n,radius):
    
    gr=nx.random_geometric_graph(n, radius, dim=2, pos=None);
    print ("Graph type: Random-Geometric")
    print ("\n number of nodes:", nx.number_of_nodes(gr))
    print ("\n number of edges:", nx.number_of_edges(gr))
    
    return gr

def generate_ws(n,k,p):
    
     return nx.connected_watts_strogatz_graph(n, k, p, tries=100, seed=None)
 
def generate_gnm(n,m):
         
    return nx.gnm_random_graph(n,m);

    
def check_bipartite(G):
    return nx.is_bipartite(G)    
 
def create_graph_from_matrix(A):

    return nx.from_numpy_matrix(A, create_using=None)    


def erdos_method(Graph,param):
    """ Initialization """
   
    adj_matrix=adjacency_matrix(Graph);  
    node_list=list_of_nodes(Graph)
  
    """ 
      1. Shuffle the node list randomly
      2. Call the function make_partition
    
    """
    
    if param == 0:
        number_of_partition=2;
        shuffled_node_list=make_shuffle(node_list)
        
        node_partition=make_partition(shuffled_node_list,number_of_partition);
    else:
        node_partition=partition_based_on_eigen_vectors(Graph,param); 
    move_flag=0;
    no_move_flag=1;
    
    
    track_movement_list=global_vertex_list(Graph);
    choice_of_partition=random.sample([0,1],1)[0];
    
    
    while no_move_flag <= 2:
    
        edge_w=compute_edges_within_partition(node_partition[choice_of_partition],adj_matrix)
        edge_a=compute_edges_across_partition(node_partition[choice_of_partition],adj_matrix)
        
    
        for vertex_u in node_partition[choice_of_partition]:
            indexvertex_u = node_partition[choice_of_partition].index(vertex_u);
            if (float(edge_w[indexvertex_u]) > float(edge_a[indexvertex_u])) and (int(track_movement_list[vertex_u])!=1):
                #print("\nMove vertex:",vertex_u, "from partition:",choice_of_partition, "to partition:", 1-int(choice_of_partition))
                move_flag=1
                node_partition=update_partition(vertex_u,node_partition,choice_of_partition);
                track_movement_list=update_global_list(vertex_u,track_movement_list);
                break;
        if move_flag==0:
           
           no_move_flag=no_move_flag + 1;
           
        else:
            
           move_flag=0;  
           no_move_flag=1;
           
        choice_of_partition=1-choice_of_partition
        
    
    ratio=float(np.sum(compute_edges_across_partition(node_partition[0],adj_matrix)))/float(Graph.number_of_edges()) 
    
    if ratio < 0.5:
      print ("something is wrong with the Erdös code..!")
      exit;
    
    return (ratio,node_partition)
    
def make_eigen_vector_based_partition(G,vec):
    
    nl= list_of_nodes(G); 
    partition_first=[];
    partition_second=[];   
    partition_final=[];
    for i in range(len(vec)):
        
        if vec[i]>0.0:
           
           partition_first.append(nl[i]);  
        
        elif vec[i]<0.0: 
           partition_second.append(nl[i]);  
        
        else:
           choice_of_partition=random.sample([0,1],1)[0];
           
           if choice_of_partition==0:
              partition_first.append(nl[i]);
           else:
              partition_second.append(nl[i]); 
               
    partition_final.append(partition_first);
    partition_final.append(partition_second);
    
    return partition_final      
            

def partition_based_on_eigen_vectors(G,mat_choice):
    
    if mat_choice==1:
        
       eig_val, eig_vec =  np.linalg.eig(adjacency_matrix(G))
       
       partition=make_eigen_vector_based_partition(G,eig_vec[:,np.nanargmin(eig_val)])
    
    if mat_choice==2:
        
       eig_val, eig_vec =  np.linalg.eig(signless_laplacian_matrix(G))
       
       partition=make_eigen_vector_based_partition(G,eig_vec[:,np.nanargmin(eig_val)])
    
    if mat_choice==3:
        
       eig_val, eig_vec =  np.linalg.eig(laplacian_matrix(G))
       
       partition=make_eigen_vector_based_partition(G,eig_vec[:,np.nanargmax(eig_val)])
           
    if mat_choice==4:
        
       eig_val, eig_vec =  np.linalg.eig(normalized_laplacian_matrix(G))
       
       partition=make_eigen_vector_based_partition(G,eig_vec[:,np.nanargmax(eig_val)])
    
    return partition   


def estrada_bipartivity(Graph):
    adj_matrix=adjacency_matrix(Graph)
    total_edge = nx.number_of_edges(Graph) 
    mine=1000;
    while not nx.is_bipartite(Graph):
        row, col=adj_matrix.shape
        bg=calculate_betaG(adj_matrix)
        for i in range(row):
            j=0;
            while j<=i:
                if adj_matrix[i,j]==1:
                   bge=calculate_betaGe(Graph,i,j)
                   e=1-(bge-bg);
                   if e<mine:
                      mine=e;
                      edgei=i;
                      edgej=j;
                j=j+1; 
        
        
        adj_matrix[edgei,edgej]=0
        adj_matrix[edgej,edgei]=0       
        Graph=nx.from_numpy_matrix(adj_matrix) 
        mine=1000;
    
    return float(nx.number_of_edges(Graph))/(total_edge)       
    
        
def calculate_betaGe(Graph,i,j):
    A=adjacency_matrix(Graph);
    A[i,j]=0;
    A[j,i]=0;
    e,v=np.linalg.eig(A)
    return sum(np.exp(-e))/sum(np.exp(e));
    
def calculate_betaG(A):
    e,v=np.linalg.eig(A)
    return sum(np.exp(-e))/sum(np.exp(e));
        
def edge_bipartivity_measure_new(Graph,param):
    
    adj_matrix=adjacency_matrix(Graph)
    total_edge = nx.number_of_edges(Graph) 
    mine=-1000;
    while not nx.is_bipartite(Graph):
        row, col=adj_matrix.shape
        for i in range(row):
            j=0;
            while j<=i:
                if adj_matrix[i,j]==1:
                   if param==1: 
                       e = edge_bipartivity_adjacency(Graph,i,j)
                   else:
                       e = edge_bipartivity_normalized_laplacian(Graph,i,j)
                   if e>mine:
                      mine=e;
                      edgei=i;
                      edgej=j;
                j=j+1; 
        
        
        adj_matrix[edgei,edgej]=0
        adj_matrix[edgej,edgei]=0       
        Graph=nx.from_numpy_matrix(adj_matrix) 
        mine=-1000;
        
    return float(nx.number_of_edges(Graph))/(total_edge) 
    

def edge_bipartivity_adjacency(Graph,i,j):
    eig_val, eig_vec =  np.linalg.eig(adjacency_matrix(Graph))
    vecmin=eig_vec[:,np.nanargmin(eig_val)]
    vecmax=eig_vec[:,np.nanargmax(eig_val)]
    
    return (vecmin[i,0]*vecmin[j,0])/(vecmax[i,0]*vecmax[j,0] + abs(vecmin[i,0]*vecmin[j,0]))

def edge_bipartivity_normalized_laplacian(Graph,i,j):    
    
    eig_val, eig_vec =  np.linalg.eig(normalized_laplacian_matrix(Graph))
    vecmax=eig_vec[:,np.nanargmax(eig_val)]
    
    return vecmax[i,0]*vecmax[j,0]
    

def read_pajek_data(pfilename):
    
    G=nx.read_pajek(pfilename);
    return G




def routine_obtain_bipartite_subgraphs(instances_of_graphs,param):
    ratio_list=[];
    for G in instances_of_graphs:
        
        is_connected_graph=check_connected(G);
         
        if is_connected_graph==True:
            """Call the function erdos_method """
            if param==0:
                tempratio_list=[];
                number_of_random_permutation=100;
                for perm in range(number_of_random_permutation):
                     [ratio,Finalnode_partition] = erdos_method(G,param) #f2
                     tempratio_list.append(ratio)
                ratio_list.append(np.max(tempratio_list))
                
            elif param==5:
                 ratio=estrada_bipartivity(G)
                 ratio_list.append(ratio)
            elif param==6:
                 ratio=edge_bipartivity_measure_new(G,1)
                 ratio_list.append(ratio) 
            elif param==7:
                 ratio=edge_bipartivity_measure_new(G,2)
                 ratio_list.append(ratio)     
            else:    
                 [ratio,Finalnode_partition] = erdos_method(G,param)
                 ratio_list.append(ratio)
        else:
             
             print("\nNot connected...")
    return ratio_list 

def call_methods(instances_of_graphs):
    ratio_erdos=routine_obtain_bipartite_subgraphs(instances_of_graphs,0)
    ratioA=routine_obtain_bipartite_subgraphs(instances_of_graphs,1)
    ratioQ=routine_obtain_bipartite_subgraphs(instances_of_graphs,2)
    ratioL=routine_obtain_bipartite_subgraphs(instances_of_graphs,3)
    ratioNL=routine_obtain_bipartite_subgraphs(instances_of_graphs,4)
    ratioEstrada=routine_obtain_bipartite_subgraphs(instances_of_graphs,5)
    ratioAnew=routine_obtain_bipartite_subgraphs(instances_of_graphs,6)
    ratioNLnew=routine_obtain_bipartite_subgraphs(instances_of_graphs,7) 
      
    return ratio_erdos,ratioA,ratioQ,ratioL, ratioNL, ratioEstrada, ratioAnew, ratioNLnew    

def save_results(erdos, rA, rQ, rL, rNL, res, rAnew, rNLnew,graph_name, outfolder):
    
    path = './'+outfolder+'/'+graph_name;
    
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    
    
    np.savetxt(path+'/ratioErdos.txt',erdos, fmt='%.3f');
    np.savetxt(path+'/ratioA.txt',rA, fmt='%.3f');
    np.savetxt(path+'/ratioQ.txt',rQ, fmt='%.3f');
    np.savetxt(path+'/ratioL.txt',rL, fmt='%.3f');
    np.savetxt(path+'/ratioNL.txt',rNL, fmt='%.3f');
    np.savetxt(path+'/ratioEstrada.txt',res, fmt='%.3f');
    np.savetxt(path+'/ratioAnew.txt',rAnew, fmt='%.3f');
    np.savetxt(path+'/ratioNLnew.txt',rNLnew, fmt='%.3f');






def histogram_plot(ratio_list,fcolor, alph,lgd):
    
    """ Controling figure properties """
    
    SMALL_SIZE = 16
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('axes', linewidth=2)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('savefig.transparent',)    
    
    #weights = np.ones_like(ratio_list)/float(len(ratio_list))
    plt.hist(ratio_list,bins=100,normed=True,facecolor=fcolor, alpha=0.4) 
    results, edges = np.histogram(ratio_list,bins=10, normed=True)
    binWidth = edges[1] - edges[0]
    plt.bar(edges[:-1], results*binWidth, binWidth,facecolor=fcolor, alpha=alph, label=lgd) 
    plt.xlim((0, 1)) 
    plt.axvline(0.5,color='r', linestyle='--',linewidth=2.0)
    plt.legend()
    plt.xlabel("ratio")
   
    
def main():
        parser = argparse.ArgumentParser()
        
        
        parser.add_argument('--g', '--graph', nargs='?', type=str,
                            default='ER', help='graph models: ER, BA, RG, WS')
        
        parser.add_argument('--p_er', '--param_er', nargs=2, type=float,
                            default=[0.0, 1.0], help='parameter for ER graph model: min max')
          
        parser.add_argument('--p_ba', '--param_ba', nargs=2, type=int,                           
                            default=[1, 10], help='parameter for BA graph model: min max')
        
        parser.add_argument('--p_rg', '--param_rg',nargs=2, type=float,                            
                            default=[0.5, 1.0], help='parameter for RG graph model: min max') 
       
        parser.add_argument('--p_ws', '--param_ws',nargs=3, type=float,                           
                            default=[8, 0.0, 0.3], help='parameters k (nearest neighbors)  &  r (probability of rewiring)for WS graph model') 
        
        parser.add_argument('--n', '--nodes',nargs=1, type=int,                           
                            default=20, help='number of nodes') 
        
        parser.add_argument('--ng', '--nog',nargs=1, type=int,                           
                            default=1, help='number of graphs') 
        
        parser.add_argument('--o', '--output',nargs=1, type=str,                           
                            default='output_folder', help='name of the output folder') 
         
        
        
        args = parser.parse_args()
        
        # Extract input parameters
        graph_name = args.g; 
        
        param_er   = args.p_er; 
        param_ba   = args.p_ba; 
        param_rg   = args.p_rg;
        param_ws   = args.p_ws;
        param_n    = args.n;
        param_ng   = args.ng;
        outfolder  = args.o;
        
        return graph_name, param_er, param_ba, param_rg, param_ws, param_n, param_ng, outfolder



if __name__ == "__main__":
    if len(sys.argv[0:])==1:
        print ("printing the default values:\n", main())
        print ("filename -h for more information\n")
    else:    
        graph_name, param_er, param_ba, param_rg, param_ws, num_nodes,number_of_graphs, outfolder=main()
        instances_of_graphs=[];
        num_graphs  =  number_of_graphs[0];
        n_nodes     = int(num_nodes[0]);
        if str(graph_name)=='ER':
           
           param_sample=np.matrix(np.random.uniform(param_er[0],param_er[1],(1,1,int(num_graphs)))) 
           for sample in range(int(num_graphs)):
               instances_of_graphs.append(generate_er(n_nodes,param_sample[0,sample]));   
        
        elif str(graph_name)=='BA':
             param_sample=np.matrix(np.random.uniform(param_ba[0],param_ba[1],(1,1,int(num_graphs)))) 
           
             for sample in range(int(num_graphs)):
                 instances_of_graphs.append(generate_ba(n_nodes,param_sample[0,sample]));

        elif str(graph_name)=='RG':
             param_sample=np.matrix(np.random.uniform(param_rg[0],param_rg[1],(1,1,int(num_graphs)))) 
           
             for sample in range(int(num_graphs)):
                 instances_of_graphs.append(generate_rg(n_nodes,param_sample[0,sample]));
 
        else:
            
             param_sample=np.matrix(np.random.uniform(param_ws[1],param_ws[2],(1,1,int(num_graphs)))) 
           
             for sample in range(int(num_graphs)):
                 instances_of_graphs.append(generate_ba(n_nodes, param_ws[0],param_sample[0,sample]));



        erdos, rA, rQ, rL, rNL, res, rAnew, rNLnew=call_methods(instances_of_graphs);  
        save_results(erdos, rA, rQ, rL, rNL, res, rAnew, rNLnew, graph_name, outfolder[0]);
