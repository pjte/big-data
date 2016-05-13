library(ggplot2)

setwd("C:/Users/Adam/Documents/GitHub/big-data/source")
library(igraph)

DE.genes <- read.csv("C:/Users/Adam/Documents/GitHub/big-data/source/shared_DEgenes.csv")
rownames(DE.genes) <- DE.genes[,2]
DE.genes <- DE.genes[,-c(1:8,68:71)]

prcomp_counts <- prcomp(t(DE.genes)) #gene wise
scores <- as.data.frame(prcomp_counts$rotation)[,c(1,2)]

set.seed(25) #make this repeatable as kmeans has random starting positions
fit <- kmeans(DE.genes, 9)
clus <- as.data.frame(fit$cluster)
head(clus)
names(clus) <- paste("cluster")

plotting <- merge(clus, scores, by = "row.names")
plotting$cluster <- as.factor(plotting$cluster)

# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 4, stat = "identity") 

#make groups based on data. This depends on the specified number of groups for the kmeans function
group1 <- subset(clus, cluster == 1)
group2 <- subset(clus, cluster == 2)
group3 <- subset(clus, cluster == 3)
group6 <- subset(clus, cluster == 6)

gene_cor <- cor(t(group6))

# cor=.3
gene_adj30 <- abs(gene_cor) > 0.3
diag(gene_adj30) <- 0
gene_graph30 <- graph.adjacency(gene_adj30, mode = "undirected") #convert adjacency to graph
comps <- clusters(gene_graph30)$membership                        #define gene cluster membership
colbar <- rainbow(max(comps)+1)                                   #define colors
V(gene_graph30)$color <- colbar[comps+1]                          #assign colors to nodes
plot(gene_graph30, layout = layout.fruchterman.reingold, vertex.size = 6, vertex.label = NA)