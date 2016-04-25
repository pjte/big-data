setwd("C:/Users/micha/Google Drive/Wyoming Classes/Data Mining/Data Mining Project/Brassica_parents/Brassica_parents")
library(data.table)
shade <- as.data.table(read.table("GC_merged_v1.5_mapping.tsv",sep = "\t", header = T))
t.shade <- as.data.table(t(shade))
setnames(t.shade,as.character(t.shade[1,]))
t.shade <- t.shade["*",:=NULL]


