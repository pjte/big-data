## Script for creating network graphs

setwd("C:/Users/micha/Google Drive/Wyoming Classes/Data Mining/Data Mining Project/Brassica_parents/Brassica_parents")
library(igraph)

shade <- read.table("GC_merged_v1.5_mapping.tsv",sep = "\t", header = T)
rownames(shade) <- shade$gene
shade <- shade[,-1]
subset.shade <- shade[sample(1:length(shade[,1]),100),]
subset.shade[is.na(subset.shade)] <- 0
gene_cor <- cor(t(subset.shade))

# cor=.3
gene_adj30 <- abs(gene_cor) > 0.3
diag(gene_adj30) <- 0
gene_graph30 <- graph.adjacency(gene_adj30, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph30)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph30)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph30, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .5
gene_adj50 <- abs(gene_cor) > 0.5
diag(gene_adj50) <- 0
gene_graph50 <- graph.adjacency(gene_adj50, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph50)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph50)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph50, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .85
gene_adj85 <- abs(gene_cor) > 0.85
diag(gene_adj85) <- 0
gene_graph85 <- graph.adjacency(gene_adj85, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph85)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph85)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph85, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

#shortest path cor=.5
gene_graph50 <- graph.adjacency(gene_adj50, mode = "undirected")
distMatrix50 <- shortest.paths(gene_graph50, v = V(gene_graph50), to = V(gene_graph50))
head(distMatrix50)[,1:7]

pl50 <- get.shortest.paths(gene_graph50, 2, 7)$vpath[[1]] # pull paths between node 2 and 7

V(gene_graph50)[pl50]$color <- paste("green")          # define node color
E(gene_graph50)$color <- paste("grey")               # define default edge color
E(gene_graph50, path = pl50)$color <- paste("blue")    # define edge color
E(gene_graph50, path = pl50)$width <- 10               # define edge width
plot(gene_graph50, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)


