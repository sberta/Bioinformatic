install.packages("igraph")
library(igraph)
setwd("C:/Users/Alessia/Desktop/università/magistrale/bioinformatica/hw2")
dat = read.csv('lccgumat.csv', header=TRUE, row.names=1, check.names=FALSE)
m = as.matrix(dat)
# build graph from adjacency matrix
g = graph.adjacency(m)
wc <- spinglass.community(g)
modularity(wc)
plot(g, vertex.label=NA, vertex.size=8, vertex.color=membership(wc),
     layout=layout.fruchterman.reingold)
wc <- louvain.community(g)
csg = cluster_spinglass(g)

