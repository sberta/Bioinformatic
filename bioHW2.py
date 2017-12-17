# -*- coding: utf-8 -*-
#%%


import os 
os.chdir("C:/Users/Alessia/Desktop/università/magistrale/bioinformatica/hw2")
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt 

sgi = pd.DataFrame.from_csv("seedGenesINTERACTOME.csv")
I = pd.DataFrame.from_csv("intersectionINTERACTOME.csv")
U = pd.DataFrame.from_csv("unionINTERACTOME.csv")            
            
#%% esercizio 1

#creating nodes
G_I = nx.from_pandas_dataframe(I, 'Protein A gene symbol', 'Protein B gene symbol')

G_sgi = nx.from_pandas_dataframe(sgi, 'Protein A gene symbol', 'Protein B gene symbol')

G_U = nx.from_pandas_dataframe(U, 'Protein A gene symbol', 'Protein B gene symbol')

#in this way i plot it
#%%
nx.draw_networkx(G_I)
plt.show()
#%%
nx.is_connected(G_I)
#%% 1.1
#number of nodes 1.1a
print(nx.number_of_nodes(G_sgi))

print(nx.number_of_nodes(G_I))

print(nx.number_of_nodes(G_U))

#%%
#number of links 1.1b
print(nx.number_of_edges(G_sgi))

print(nx.number_of_edges(G_I))
    
print(nx.number_of_edges(G_U))

#%%
#number of isolated nodes 1.1c
nodeIsolGI = nx.isolates(G_I)

nodeIsolGU = nx.isolates(G_U)

nodeIsolSGI = nx.isolates(G_sgi)

#%%
#average path length 1.1d

#asplGI = nx.average_shortest_path_length(G_I)

asplGU = nx.average_shortest_path_length(G_U)

asplsgi = nx.average_shortest_path_length(G_sgi)

#%%
#average path length 1.1e

adcGI = nx.average_degree_connectivity(G_I)

adcGU = nx.average_degree_connectivity(G_U)

adcSGI = nx.average_degree_connectivity(G_sgi)

#%%
#average clustering 1.1f

avgclustGI = nx.average_clustering (G_I)

avgclustGU = nx.average_clustering (G_U)

avgclustSGI = nx.average_clustering (G_sgi)

#%%
#average clustering 1.1g

#diameterGI = nx.diameter(G_I)

diameterGU =nx.diameter(G_U)
radGU =nx.radius(G_U)

diameterSGI =nx.diameter(G_sgi)

#%%
import operator
centrGU = nx.degree_centrality(G_U)
print(max(centrGU, key=centrGU.get))

centrGI = nx.degree_centrality(G_I)
print(max(centrGU, key=centrGU.get))

centrSGI = nx.degree_centrality(G_sgi)
print(max(centrGU, key=centrGU.get))

#%%1.2

LCCgi = max(nx.connected_component_subgraphs(G_I), key=len)
lccgimat = nx.to_pandas_dataframe(LCCgi)
lccgimat.to_csv("lccgimat.csv")
LCCgu = max(nx.connected_component_subgraphs(G_U), key=len)
lccgumat = nx.to_pandas_dataframe(LCCgu)
lccgumat.to_csv("lccgumat.csv")
#LCCsgi = max(nx.connected_component_subgraphs(G_sgi), key=len)

#%% 1.2 a 1
#number of nodes 
print(nx.number_of_nodes(LCCgi))

print(nx.number_of_nodes(LCCgu))

#print(nx.number_of_nodes(LCCsgi))

#%%
#number of links 1.2 a 2
print(nx.number_of_edges(LCCgi))

#print(nx.number_of_edges(LCCsgi))
    
print(nx.number_of_edges(LCCgu))


#%%
#average path length 1.2 a 3

asplLCCgi = nx.average_shortest_path_length(LCCgi)

asplLCCgu = nx.average_shortest_path_length(LCCgu)

#asplLCCsgi = nx.average_shortest_path_length(LCCsgi)

#%%
#average path length 1.2 a 4

adcLCCgi = nx.average_degree_connectivity(LCCgi)

adcLCCgu = nx.average_degree_connectivity(LCCgu)

#adcLCCsgi = nx.average_degree_connectivity(LCCsgi)

#%%
#average clustering 1.2 a 4

avgclustLCCgi = nx.average_clustering (LCCgi)

avgclustLCCgu = nx.average_clustering (LCCgu)

#avgclustLCCsgi = nx.average_clustering (LCCsgi)

#%%
#average clustering 1.2 a 5

diameterLCCgi = nx.diameter(LCCgi)

diameterLCCgu =nx.diameter(LCCgu)

#diameterLCCsgi =nx.diameter(LCCsgi)

#%% 1.2 a 6
import operator
centrLCCgi = nx.degree_centrality(LCCgi)
print(max(centrLCCgi, key=centrLCCgi.get))

centrLCCgu = nx.degree_centrality(LCCgu)
print(max(centrLCCgu, key=centrLCCgu.get))
#
#centrLCCsgi = nx.degree_centrality(LCCsgi)
#print(max(centrLCCsgi, key=centrLCCsgi.get))


#%%1.2 b1

degreeLCCgi = nx.degree(LCCgi)


degreeLCCgu = nx.degree(LCCgu)

#
#degreeLCCsgi = nx.degree(LCCsgi)


#%%1.2 b2

betwLCCgi = nx.betweenness_centrality(LCCgi, normalized=True)

#betwLCCsgi = nx.betweenness_centrality(LCCsgi, normalized=True)

betwLCCgu = nx.betweenness_centrality(LCCgu, normalized=True)

#%% 1.2 b3

eig_centrLCCgi = nx.eigenvector_centrality_numpy(LCCgi) #dava errore
#%%
#eig_centrLCCsgi = nx.eigenvector_centrality(LCCsgi)

eig_centrLCCgu = nx.eigenvector_centrality(LCCgu)
#%% 1.2 b4

clos_cenLCCgi = nx.closeness_centrality(LCCgi, normalized=True)

#clos_cenLCCsgi = nx.closeness_centrality(LCCsgi, normalized=True)

clos_cenLCCgu = nx.closeness_centrality(LCCgu, normalized=True)

#%% 1.2 b5    


nodeRatioLCCgi =  {k: betwLCCgi[k]/degreeLCCgi[k] for k in betwLCCgi.keys()} 

#nodeRatioLCCsgi = {k: betwLCCsgi[k]/degreeLCCsgi[k] for k in betwLCCsgi.keys()}  

nodeRatioLCCgu = {k: betwLCCgu[k]/degreeLCCgu[k] for k in betwLCCgu.keys()}               
                       


#%% 2


import markov_clustering as mc
import networkx as nx
import random
 #%%
# number of nodes to use
numnodes = nx.number_of_nodes(LCCgi)

# then get the adjacency matrix (in sparse form)
matrix = nx.to_scipy_sparse_matrix(LCCgi)
#%%
result = mc.run_mcl(matrix)            # run MCL with default parameters
clusters = mc.get_clusters(result)     # get clusters

mc.draw_graph(matrix, clusters,  node_size=300, with_labels=False, edge_color="silver")
#%%


 #%%
# number of nodes to use
numnodes = nx.number_of_nodes(LCCgu)

# then get the adjacency matrix (in sparse form)
matrix = nx.to_scipy_sparse_matrix(LCCgu)
#%%
result = mc.run_mcl(matrix)            # run MCL with default parameters
clusters = mc.get_clusters(result)     # get clusters

mc.draw_graph(matrix, clusters,  node_size=500, with_labels=False, edge_color="silver")
  
#%%
#Louvain - This requires python-louvain library do "pip install python-louvain"
import community
community.best_partition(LCCgu)

#ploting louvain
partition = community.best_partition(LCCgu)

size = float(len(set(partition.values())))
pos = nx.spring_layout(LCCgu)
count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(LCCgu, pos, list_nodes, node_size = 20,
                                node_color = str(count / size))


nx.draw_networkx_edges(LCCgi,pos, alpha=0.5)
plt.show()             

#%%

import community
import networkx as nx
import matplotlib.pyplot as plt

#better with karate_graph() as defined in networkx example.
#erdos renyi don't have true community structure


#first compute the best partition
partition = community.best_partition(LCCgu)
values = [partition.get(node) for node in LCCgu.nodes()]

nx.draw_spring(LCCgu,cmap = plt.get_cmap('jet'), node_color = values, node_size=700)

#%%
plt.figure()
partition = community.best_partition(LCCgi)
values = [partition.get(node) for node in LCCgi.nodes()]

nx.draw_spring(LCCgi,cmap = plt.get_cmap('jet'), node_color = values, node_size=700)

#%%

wc = spinglass.community(LCCgu)

