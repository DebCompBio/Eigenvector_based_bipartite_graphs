#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 17:29:24 2018

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
    
    NOTE: the ratio should always be greater than or equal to 0.5 
"""


""" import packages """

import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import time


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

def  LaplacianMatrix(G):
    
      return np.matrix(DegreeMatrix(G) - AdjacencyMatrix(G)) 

def SignLessLaplacianMatrix(G):

      return np.matrix(DegreeMatrix(G) + AdjacencyMatrix(G))

def NormaLizedLaplacianMatrix(G):
    
     DiagSqrt = DegreeMatrix(G);
     
     for row in range(len(DiagSqrt)):
         
         DiagSqrt[row,row] = pow(DiagSqrt[row,row],-0.5);
         
     return np.matrix(DiagSqrt*LaplacianMatrix(G)*DiagSqrt)    
 
def DegreeMatrix(G):
     A=AdjacencyMatrix(G) 
     trA=np.transpose(np.sum(A,1))
     dye=np.eye(len(G.nodes()));
     for dd in range(len(G.nodes())):
         
         dye[dd,dd]=trA[0,dd];
     
     return dye
         

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


def CheckConnected(G):
    
    """Check if the graph is connected or not """
    
    return nx.is_connected(G);
    
def EdwardsRatio(G):
    
   return float(G.number_of_edges())/float(2) + (float(len(G.nodes()))-1.0)/4.0;
    

""" Graph generation routines """

def GenerateER(n,p):
    
    return nx.erdos_renyi_graph(n, p, seed=None, directed=False);

def GenerateBA(n,m):
    
    return  nx.barabasi_albert_graph(n, m, seed=None);
    
def GenerateRG(n,radius):
    
    return nx.random_geometric_graph(n, radius, dim=2, pos=None);

def GenerateWS(n,k,p):
    
     return nx.connected_watts_strogatz_graph(n, k, p, tries=100, seed=None)
 
    
def CheckBipartite(G):
    return nx.is_bipartite(G)    
 
def CreateGraphFromMatrix(A):

    return nx.from_numpy_matrix(A, create_using=None)    


def ErdosMethod(Graph,param):
    """ Initialization """
   
    AdjMatrix=AdjacencyMatrix(Graph);  
    NodeList=ListOfNodes(Graph)
  
    """ 
      1. Shuffle the node list randomly
      2. Call the function MakePartition
    
    """
    
    if param == 0:
        NumberOfPartition=2;
        ShuffledNodeList=MakeShuffle(NodeList)
        
        NodePartition=MakePartition(ShuffledNodeList,NumberOfPartition);
    else:
        NodePartition=EigenVectors(Graph,param); 
    move_flag=0;
    no_move_flag=1;
    
    # Create a global list to track the movement of vertices
    
    TrackMovementList=GlobalVertexList(Graph);
    ChoiceOfPartition=random.sample([0,1],1)[0];
    
    
    while no_move_flag <= 2:
    
        EdgeW=ComputeEdgesWithinPartition(NodePartition[ChoiceOfPartition],AdjMatrix)
        EdgeA=ComputeEdgesAcrossPartition(NodePartition[ChoiceOfPartition],AdjMatrix)
        
    
        for vertexU in NodePartition[ChoiceOfPartition]:
            indexVertexU = NodePartition[ChoiceOfPartition].index(vertexU);
            if (float(EdgeW[indexVertexU]) > float(EdgeA[indexVertexU])) and (int(TrackMovementList[vertexU])!=1):
                """ Move a vertex """
                #print("\nMove vertex:",vertexU, "from partition:",ChoiceOfPartition, "to partition:", 1-int(ChoiceOfPartition))
                move_flag=1
                NodePartition=UpdatePartition(vertexU,NodePartition,ChoiceOfPartition);
                TrackMovementList=UpdateGlobalList(vertexU,TrackMovementList);
                break;
        if move_flag==0:
           
           no_move_flag=no_move_flag + 1;
           
        else:
            
           move_flag=0;  
           no_move_flag=1;
        """ Change the starting partition """
        ChoiceOfPartition=1-ChoiceOfPartition
        
    
    ratio=float(np.sum(ComputeEdgesAcrossPartition(NodePartition[0],AdjMatrix)))/float(Graph.number_of_edges()) 
    
    if ratio < 0.5:
      print ("something is wrong with the Erdös code..!")
      exit;
    
    return (ratio,NodePartition)
    
def MakeEigenVectorBasedPartition(G,vec):
    
    nl= ListOfNodes(G); 
    part1=[];
    part2=[];   
    FinalPart=[];
    for i in range(len(vec)):
        
        if vec[i]>0.0:
           
           part1.append(nl[i]);  
        
        elif vec[i]<0.0: 
           part2.append(nl[i]);  
        
        else:
           ChoiceOfPartition=random.sample([0,1],1)[0];
           
           if ChoiceOfPartition==0:
              part1.append(nl[i]);
           else:
              part2.append(nl[i]); 
               
    FinalPart.append(part1);
    FinalPart.append(part2);
    
    return FinalPart      
            

def EigenVectors(G,mat_choice):
    
    if mat_choice==1:
        
       EigVal, EigVec =  np.linalg.eig(AdjacencyMatrix(G))
       
       Partition=MakeEigenVectorBasedPartition(G,EigVec[:,np.nanargmin(EigVal)])
    
    if mat_choice==2:
        
       EigVal, EigVec =  np.linalg.eig(SignLessLaplacianMatrix(G))
       
       Partition=MakeEigenVectorBasedPartition(G,EigVec[:,np.nanargmin(EigVal)])
    
    if mat_choice==3:
        
       EigVal, EigVec =  np.linalg.eig(LaplacianMatrix(G))
       
       Partition=MakeEigenVectorBasedPartition(G,EigVec[:,np.nanargmax(EigVal)])
           
    if mat_choice==4:
        
       EigVal, EigVec =  np.linalg.eig(NormaLizedLaplacianMatrix(G))
       
       Partition=MakeEigenVectorBasedPartition(G,EigVec[:,np.nanargmax(EigVal)])
    
    return Partition   


def EstradaBipartivity(Graph):
    AdjMatrix=AdjacencyMatrix(Graph)
    TotalEdge = nx.number_of_edges(Graph) 
    mine=1000;
    while not nx.is_bipartite(Graph):
        row, col=AdjMatrix.shape
        bg=CalculateBetaG(AdjMatrix)
        for i in range(row):
            j=0;
            while j<=i:
                if AdjMatrix[i,j]==1:
                   bge=CalculateBetaGe(Graph,i,j)
                   e=1-(bge-bg);
                   if e<mine:
                      mine=e;
                      edgei=i;
                      edgej=j;
                j=j+1; 
        
        
        AdjMatrix[edgei,edgej]=0
        AdjMatrix[edgej,edgei]=0       
        Graph=nx.from_numpy_matrix(AdjMatrix) 
        mine=1000;
    
    return float(nx.number_of_edges(Graph))/(TotalEdge)       
    
        
def CalculateBetaGe(Graph,i,j):
    A=AdjacencyMatrix(Graph);
    A[i,j]=0;
    A[j,i]=0;
    e,v=np.linalg.eig(A)
    return sum(np.exp(-e))/sum(np.exp(e));
    
def CalculateBetaG(A):
    e,v=np.linalg.eig(A)
    return sum(np.exp(-e))/sum(np.exp(e));
        
def NewEdgeBipartivitMeasure(Graph,param):
    
    AdjMatrix=AdjacencyMatrix(Graph)
    TotalEdge = nx.number_of_edges(Graph) 
    mine=-1000;
    while not nx.is_bipartite(Graph):
        row, col=AdjMatrix.shape
        for i in range(row):
            j=0;
            while j<=i:
                if AdjMatrix[i,j]==1:
                   if param==1: 
                       e = AMatrixBasedMeasure(Graph,i,j)
                   else:
                       e = NLMatrixBasedMeasure(Graph,i,j)
                   if e>mine:
                      mine=e;
                      edgei=i;
                      edgej=j;
                j=j+1; 
        
        
        AdjMatrix[edgei,edgej]=0
        AdjMatrix[edgej,edgei]=0       
        Graph=nx.from_numpy_matrix(AdjMatrix) 
        mine=-1000;
        
    return float(nx.number_of_edges(Graph))/(TotalEdge) 
    

def AMatrixBasedMeasure(Graph,i,j):
    EigVal, EigVec =  np.linalg.eig(AdjacencyMatrix(Graph))
    vecmin=EigVec[:,np.nanargmin(EigVal)]
    vecmax=EigVec[:,np.nanargmax(EigVal)]
    
    return (vecmin[i,0]*vecmin[j,0])/(vecmax[i,0]*vecmax[j,0] + abs(vecmin[i,0]*vecmin[j,0]))

def NLMatrixBasedMeasure(Graph,i,j):    
    
    EigVal, EigVec =  np.linalg.eig(NormaLizedLaplacianMatrix(Graph))
    vecmax=EigVec[:,np.nanargmax(EigVal)]
    
    return vecmax[i,0]*vecmax[j,0]
    

    

""" 
    Define a graph 

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
def ReadPajekData(pfilename):
    
    G=nx.read_pajek(pfilename);
    return G




def MainRoutine(GraphInstance,param):
    ratioList=[];
    for G in GraphInstance:
        
        IsConnectedGraph=CheckConnected(G);
         
        if IsConnectedGraph==True:
            """Call the function ErdosMethod """
            if param==0:
                tempratioList=[];
                NumberOfRandomPermutation=100;
                for perm in range(NumberOfRandomPermutation):
                     [ratio,FinalNodePartition] = ErdosMethod(G,param) #f2
                     tempratioList.append(ratio)
                ratioList.append(np.max(tempratioList))
                
            elif param==5:
                 ratio=EstradaBipartivity(G)
                 ratioList.append(ratio)
            elif param==6:
                 ratio=NewEdgeBipartivitMeasure(G,1)
                 ratioList.append(ratio) 
            elif param==7:
                 ratio=NewEdgeBipartivitMeasure(G,2)
                 ratioList.append(ratio)     
            else:    
                 [ratio,FinalNodePartition] = ErdosMethod(G,param)
                 ratioList.append(ratio)
        else:
             
             print("\nNot connected...")
    return ratioList 

def PlotHistogram(ratioList,fcolor, alph,lgd):
    
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
    
    #weights = np.ones_like(ratioList)/float(len(ratioList))
    plt.hist(ratioList,bins=100,normed=True,facecolor=fcolor, alpha=0.4) 
    results, edges = np.histogram(ratioList,bins=10, normed=True)
    binWidth = edges[1] - edges[0]
    plt.bar(edges[:-1], results*binWidth, binWidth,facecolor=fcolor, alpha=alph, label=lgd) 
    plt.xlim((0, 1)) 
    plt.axvline(0.5,color='r', linestyle='--',linewidth=2.0)
    plt.legend()
    plt.xlabel("ratio")
   
    
NumberOfGraphs=1200;    

#a=np.matrix(np.random.uniform(0.2,1,(1,1,NumberOfGraphs)))   #ER model param

#m = np.matrix(np.random.randint(1,10,(1,1,NumberOfGraphs)))  #BA model param

#r = np.matrix(np.random.uniform(0.5,1,(1,1,NumberOfGraphs))) #RG model threshold

prew = np.matrix(np.random.uniform(0,0.3,(1,1,NumberOfGraphs))) # WS -model rewiring probability

ratioList=[];
start_time = time.clock()
GraphInstance=[];
GraphModelName='WSmodel/4' 
   
for pr in range(NumberOfGraphs):
        
        #GraphInstance.append(GenerateER(20,a[0,pr]));     #ER Graphmodel
        
        #GraphInstance.append(GenerateBA(20,m[0,pr]));     #BA Graphmodel
        
        #GraphInstance.append(GenerateRG(20,r[0,pr]));     #RG Graphmodel
        
        GraphInstance.append(GenerateWS(20,4,prew[0,pr]));        #WS Graphmodel
        
################### Erdös method ###################
print("\nEntering Erdos method")
ratioerdos=MainRoutine(GraphInstance,0)
np.savetxt('./'+GraphModelName+'/ratioErdos.txt',ratioerdos, fmt='%.3f'); 
#PlotHistogram(ratioerdos,'r',0.2,'Erdös method'); 

print("\nEntering Eigenvector based methods")
################### Adjacency matrix + Erdös method ###################
ratioA=MainRoutine(GraphInstance,1)
#PlotHistogram(ratioA,'g',0.2,'A-matrix'); 
np.savetxt('./'+GraphModelName+'/ratioA.txt',ratioA, fmt='%.3f'); 

################### Signless Laplacian matrix + Erdös method ###################
ratioQ=MainRoutine(GraphInstance,2)
#PlotHistogram(ratioQ,'b',0.2,'Q-matrix'); 
np.savetxt('./'+GraphModelName+'/ratioQ.txt',ratioQ, fmt='%.3f'); 

################### Laplacian matrix + Erdös method ###################
ratioL=MainRoutine(GraphInstance,3)
#PlotHistogram(ratioL,'k',0.2,'Q-matrix'); 
np.savetxt('./'+GraphModelName+'/ratioL.txt',ratioL, fmt='%.3f'); 

################### normalized Laplacian matrix + Erdös method ###################
ratioNL=MainRoutine(GraphInstance,4)
#PlotHistogram(ratioNL,'y',0.2,'NL-matrix'); 
np.savetxt('./'+GraphModelName+'/ratioNL.txt',ratioNL, fmt='%.3f'); 

print("\nEntering Estrada edge bipartivity")
################### Estrada edge bipartivity ###################
ratioEstrada=MainRoutine(GraphInstance,5)
#PlotHistogram(ratioEstrada,'y',0.2,'Estrada-edge'); 
np.savetxt('./'+GraphModelName+'/ratioEstrada.txt',ratioEstrada, fmt='%.3f'); 

print("\nEntering New measures")
################### New measure based on the Adjacency matrix ###################
ratioAnew=MainRoutine(GraphInstance,6)
#PlotHistogram(ratioAnew,'y',0.2,'A-new'); 
np.savetxt('./'+GraphModelName+'/ratioAnew.txt',ratioAnew, fmt='%.3f'); 

################### New measure based on the Normalized Laplacian matrix ###################
ratioNLnew=MainRoutine(GraphInstance,7)
#PlotHistogram(ratioNLnew,'y',0.2,'NL-new'); 
np.savetxt('./'+GraphModelName+'/ratioNLnew.txt',ratioNLnew, fmt='%.3f'); 

plt.show()
 
print (time.clock() - start_time, "seconds")