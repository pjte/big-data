## Script for creating network graphs
setwd("C:/Users/micha/Google Drive/Wyoming Classes/Data Mining/Data Mining Project/Brassica_parents/Brassica_parents")
library(igraph)

# shade <- read.table("GC_merged_v1.5_mapping.tsv",sep = "\t", header = T)
DE.genes <- read.table("C:/Users/micha/Google Drive/Wyoming Classes/Data Mining/Data Mining Project/big-data/source/shared_DEgenes.csv",sep = ",", header = T)
# DE.genes.names <- DE.genes$V1
# GxE_counts <- as.data.frame(shade[DE.genes.names,])
rownames(DE.genes) <- DE.genes[,2]
DE.genes <- DE.genes[,-c(1:8,68:71)]
genes_cor <- cor(t(DE.genes)) # calculate the correlation between all gene pairs

# cor=.3
gene_adj30 <- abs(genes_cor) > 0.3
diag(gene_adj30) <- 0
gene_graph30 <- graph.adjacency(gene_adj30, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph30)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph30)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph30, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .5
gene_adj50 <- abs(genes_cor) > 0.5
diag(gene_adj50) <- 0
gene_graph50 <- graph.adjacency(gene_adj50, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph50)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph50)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph50, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .70
gene_adj70 <- abs(genes_cor) > 0.70
diag(gene_adj70) <- 0
gene_graph70 <- graph.adjacency(gene_adj70, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph70)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph70)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph70, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .80
gene_adj80 <- abs(genes_cor) > 0.80
diag(gene_adj80) <- 0
gene_graph80 <- graph.adjacency(gene_adj80, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph80)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph80)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph80, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .85
gene_adj85 <- abs(genes_cor) > 0.85
diag(gene_adj85) <- 0
gene_graph85 <- graph.adjacency(gene_adj85, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph85)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph85)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph85, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .90
gene_adj90 <- abs(genes_cor) > 0.90
diag(gene_adj90) <- 0
gene_graph90 <- graph.adjacency(gene_adj90, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph90)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph90)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph90, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .95
gene_adj95 <- abs(genes_cor) > 0.95
diag(gene_adj95) <- 0
gene_graph95 <- graph.adjacency(gene_adj95, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph95)$membership                        #define gene cluster membership
group1 <- subset(comps, comps == 1)
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph95)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph95, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

# cor= .95
# only group 1
gene_adj95 <- abs(genes_cor[names(group1),names(group1)]) > 0.95
diag(gene_adj95) <- 0
gene_graph95 <- graph.adjacency(gene_adj95, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph95)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph95)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph95, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)

#shortest path cor=85
gene_graph85 <- graph.adjacency(gene_adj85, mode = "undirected")
distMatrix85 <- shortest.paths(gene_graph85, v = V(gene_graph85), to = V(gene_graph85))
head(distMatrix85)[,1:7]

pl85 <- get.shortest.paths(gene_graph85, 2, 7)$vpath[[1]] # pull paths between node 2 and 7

V(gene_graph85)[pl85]$color <- paste("green")          # define node color
E(gene_graph85)$color <- paste("grey")               # define default edge color
E(gene_graph85, path = pl85)$color <- paste("blue")    # define edge color
E(gene_graph85, path = pl85)$width <- 10               # define edge width
plot(gene_graph85, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)


