# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 22:50:45 2018

@author: Shiva
"""

# Python program to print all paths from a source to destination. 
import numpy as np
import networkx as nx

def Fault_Isolation(F):
   
    switches=np.matrix([[54, 94], [117, 123], [1, 125], [97, 124], [18, 116], [60, 119], [13, 121]]); 
    # Create a graph from the given edges
    nNodes=125;
    nEdges=126;
    edges = np.loadtxt('edges.txt'); 
    
    n=(nNodes+1,nNodes+1)   
    adjMatrix = np.zeros(n)
         
        # scan the arrays edge_u and edge_v
    for i in range(0, nEdges):
        u=int(edges[i,0])
        v=int(edges[i,1])
        adjMatrix[u,v]=1;
    
    G=nx.from_numpy_matrix(adjMatrix)
    Source=125
    Fault=F;
    
    # ways stores all possible paths from source to fault location
    ways=list(nx.all_simple_paths(G,source=Source,target=Fault))
    #print(ways) # all simple paths
    
    vector=[];
    for i in range(0,ways.__len__()):
        ls=ways[i].__len__()   
        a=ways[i]
        isolate=[]; 
        for m in range(0, ls-1):
                store=[a[ls-m-1], a[ls-m-2]]
                #print(store)
                # check if the line in store is switch. If yes, then we need to open that line for fault isoltation. 
                for k in range(0, switches.__len__()):
                    if store[0]==switches[k,0]:
                        if store[1]==switches[k,1]:
                            isolate=switches[k]                        
                            vector.append(isolate)
                            #print(isolate)
                                           
                    if store[0]==switches[k,1]:
                        if store[1]==switches[k,0]:
                            isolate=switches[k]                        
                            vector.append(isolate)
                            #print(isolate)
                if isolate.__len__()==1:
                        break
    
    open_lines=[]
    isolate_Fault=[]  
    vector=np.matrix(np.array(vector))
    print('\n \n The edgelist to open are:')
    print(vector)
    
    for k in range (0, vector.__len__()):
        inde=np.where(np.all(vector[k]==edges,axis=1))
        open_lines.append(inde)
        a=np.array(open_lines)
        isolate_Fault.append(a[k,0,0])
        #print(a[k,0,0])
    
    isolate_Fault=set(isolate_Fault)
    print('\n \n The line index for isolating the faults are:')
    print(list(isolate_Fault))
    
    open_sw=list(isolate_Fault)
    return open_sw
                    
            
