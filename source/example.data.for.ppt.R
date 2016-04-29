setwd("C:/Users/micha/Google Drive/Wyoming Classes/Data Mining/Data Mining Project/Brassica_parents/Brassica_parents")
library(data.table)
shade <- read.table("GC_merged_v1.5_mapping.tsv",sep = "\t", header = T)

example.data <- shade[c(2:13,15:31),
c(
"gene",
"IMB211_SHADE_1_LEAF.1.merged.bam",
"IMB211_SHADE_1_ROOT.1.merged.bam",
"IMB211_SUN_1_LEAF.1.merged.bam",
"IMB211_SUN_1_ROOT.1.merged.bam",
"R500_SHADE_1_LEAF.1.merged.bam",            
"R500_SHADE_1_ROOT.1.merged.bam",            
"R500_SUN_1_LEAF.1.merged.bam",              
"R500_SUN_1_ROOT.1.merged.bam"              
)]

names(example.data)<-
  c(
    "gene",
    "IMB211_SHADE_1_LEAF",
    "IMB211_SHADE_1_ROOT",
    "IMB211_SUN_1_LEAF",
    "IMB211_SUN_1_ROOT",
    "R500_SHADE_1_LEAF",            
    "R500_SHADE_1_ROOT",            
    "R500_SUN_1_LEAF",              
    "R500_SUN_1_ROOT"
  )

IMB211_leaf_diff <- example.data$IMB211_SUN_1_LEAF - example.data$IMB211_SHADE_1_LEAF
IMB211_root_diff <- example.data$IMB211_SUN_1_ROOT - example.data$IMB211_SHADE_1_ROOT
IMB211_root_diff <- example.data$R500_SUN_1_LEAF - example.data$R500_SHADE_1_LEAF
r500_root_diff <- example.data$R500_SUN_1_ROOT - example.data$R500_SHADE_1_ROOT

example.differential.data <- cbind(IMB211_leaf_diff,IMB211_root_diff,
                                   IMB211_root_diff,r500_root_diff)
example.differential.data[is.na(example.differential.data)] <- 0
row.names(example.differential.data) <- example.data$gene

View(t(example.differential.data))





