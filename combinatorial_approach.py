#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 17:29:24 2018

@author: Paul, Debdas 

"""
""" Description 
    
    The script COMBINATORIAL_APPROACH is an algorithmic implementation of Erdos'
    result on making a graph bipartite described in:
    
        Erdös P. On some extremal problems in graph theory, 
                 Israel Journal of Matematics,  (1965); 3(2):pp. 113-6
    
    INPUT :  A simple, connected, and undirected graph G
    
    OUTPUT: A bipartiiton of vertices with the ratio between the edges in the bipartition
            to the total number of edges in the graph G
    
    NOTE: the ration should always be greater than or equal to 0.5 
"""


""" import packages """

import networkx as nx
import random
import numpy as np

def MakePartition(lst,n):
    """ Make random partition of vertices """
    
    division = len(lst) / float(n) 
    return [lst[int(round(division * i)):
            int(round(division * (i + 1)))] for i in range(n)] 

def MakeShuffle(lst):  

    return random.sample(lst,len(lst));      
    
def UpdatePartition(vertexU,NodePartition,ChoiceOfPartition):

    """ Move a vertex from one partition to another"""
    NodePartition[ChoiceOfPartition].remove(vertexU)
    NodePartition[1-int(ChoiceOfPartition)].extend([vertexU])
    return NodePartition
    
 
    
def GlobalVertexList(G):
   TrackMovementList={};
   
   for node in  ListOfNodes(G):
       
       TrackMovementList[node]=0;
       
   return TrackMovementList 

def UpdateGlobalList(vertexU,TrackMovementList):
    
    """ Update the global list of vertices to check
        if there will be further movement of vertices"""
        
    TrackMovementList[vertexU]=1;
    return TrackMovementList;
       
def  AdjacencyMatrix(G):
    
    return nx.to_numpy_matrix(G);

def ListOfNodes(G):
    
    """ Return the list of nodes of graph G """
    
    return list(G.nodes)

def ComputeEdgesWithinPartition(partition,AdjMatrix):
    """ Compute the neighbors of a vertex u within its member 
        partion"""
    
    pr=AdjMatrix[np.ix_(partition, partition)] 
    
    return np.sum(pr,axis=1)
    
    
def ComputeEdgesAcrossPartition(partition,AdjMatrix):
    """ Compute the neighbors of a vertex u across its member 
        partion"""
    return ComputeEdges(AdjMatrix)[np.ix_(partition, [0])]-ComputeEdgesWithinPartition(partition,AdjMatrix);

def ComputeEdges(Adj):
    
    return np.sum(Adj,axis=1)
    

def ErdosMethod(Graph):
    """ Initialization """
    #Graph=CreateGraph()
    
    AdjMatrix=AdjacencyMatrix(Graph);
    
    print("The graph adjacency matrix is:\n",AdjMatrix)       
    NodeList=ListOfNodes(Graph)
    
    NumberOfPartition=2;

    """ 
      1. Shuffle the node list randomly
      2. Call the function MakePartition
    
    """
    
    ShuffledNodeList=MakeShuffle(NodeList)
#    
    NodePartition=MakePartition(ShuffledNodeList,NumberOfPartition);
    
    move_flag=0;
    no_move_flag=1;
    
    # Create a global list to track the movement of vertices
    
    TrackMovementList=GlobalVertexList(Graph);
    ChoiceOfPartition=random.sample([0,1],1)[0];
    while no_move_flag < 2:
        print ("\nNode partition: ", NodePartition)
    
        EdgeW=ComputeEdgesWithinPartition(NodePartition[ChoiceOfPartition],AdjMatrix)
        EdgeA=ComputeEdgesAcrossPartition(NodePartition[ChoiceOfPartition],AdjMatrix)
        
        #print(EdgeW)
        
        #print("\n",EdgeA)    
    
        for vertexU in NodePartition[ChoiceOfPartition]:
            indexVertexU = NodePartition[ChoiceOfPartition].index(vertexU);
            
            if (float(EdgeW[indexVertexU]) > float(EdgeA[indexVertexU])) and (TrackMovementList[vertexU]!=1):
                """ Move a vertex """
                print("\nMove vertex:",vertexU, "from partition:",ChoiceOfPartition, "to partition:", 1-int(ChoiceOfPartition))
                move_flag=1
                NodePartition=UpdatePartition(vertexU,NodePartition,ChoiceOfPartition);
                TrackMovementList=UpdateGlobalList(vertexU,TrackMovementList);
                #print(NodePartition)
                break;
        if move_flag==0:
           
           no_move_flag=no_move_flag + 1;
           
        else:
            
           move_flag=0;  
           
        """ Change the starting partition """
        ChoiceOfPartition=1-ChoiceOfPartition
    
    print("\nNumber of edges in the bipartition: ", np.sum(ComputeEdgesAcrossPartition(NodePartition[0],AdjMatrix)))
    
    print ("\nVertex Movement list : ",TrackMovementList)
    ratio=float(np.sum(ComputeEdgesAcrossPartition(NodePartition[0],AdjMatrix)))/float(Graph.number_of_edges()) 
    
    return (ratio,NodePartition)

""" Define a graph """

G = nx.petersen_graph(); # petersen graph, modify here according to your choice of the graph
    
"""Call the function ErdosMethod """

[ratio,FinalNodePartition] = ErdosMethod(G)

print ("\nFor petersen graph the ratio is: ", ratio, "with partition: ", FinalNodePartition)

EdwardsRatio = float(G.number_of_edges())/float(2) + (float(len(G.nodes()))-1.0)/4.0;

print ("\n Edwards ratio for the graph is: " ,EdwardsRatio/G.number_of_edges())
    