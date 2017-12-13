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
nx.draw(G_sgi)
plt.show
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

diameterSGI =nx.diameter(G_sgi)

#%%

centrGU = nx.degree_centrality(G_U)

centrGI = nx.degree_centrality(G_I)

centrSGI = nx.degree_centrality(G_sgi)






