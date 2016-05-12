source("https://bioconductor.org/biocLite.R")
biocLite("baySeq")
biocLite("edgeR")
library(edgeR)

# Read in data file
dat <- read.table(file = "/Users/mblai/Documents/BigData/Brassica_parents/GC_merged_v1.5_mapping.tsv",
                  header = T)

# Remove * row and make gene IDs row labels
data <- dat[2:41921, 2:67]
rownames(data) <- dat[2:41921,1]
# Shorten the names
colnames(data)<-sub(".1.merged.bam", "", colnames(data))

# Convert NAs to zeros
data[is.na(data)]<-0
# Round decimals up to nearest integer
data<-ceiling(data)

#remove apical meristem
data<-data[,c(-1,-7,-32,-39,-51)]
#remove floral meristem
data<-data[,c(-30,-36)]


# Filter out lowly expressed genes
#  a. check distribution of counts for unfiltered genes
summary(rowSums(data))
#  b. filter genes > 10 reads in >= 3 samples
filtered <- data[rowSums(data > 10) >= 3,]
#  c. check distribution of counts for filtered genes
summary(rowSums(filtered))

# Boxplots - Added 1 to all filtered counts to avoid log2(0)
png("Final_pre-normalized_boxplots.png")
boxplot(log2(filtered+1),ylab="log2(Counts)",
        main="Pre-Normalized Filtered RNA-Seq Counts",xlab="Sample")
# Two very interesting samples exist
dev.off()

# Scatterplot
png("Final_corr_scatterplot.png")
pairs(log2(filtered[1:1000,1:10]))
dev.off()
# Notice biological replicates have high correlation

# Assign groups based on combinations of genotype and treatment
# And put in a dataframe for edgeR
grouptype <- data.frame(
  sample=colnames(filtered),
  geno=regmatches(colnames(filtered),regexpr("R500|IMB211",colnames(filtered))),
  trmt=regmatches(colnames(filtered),regexpr("SHADE|SUN",colnames(filtered)))
)
grouptype$group <- paste(grouptype$geno,grouptype$trmt,sep="_")

# Set the reference treatment to "SUN"
grouptype$trmt <- relevel(grouptype$trmt,ref="SUN")

# Created DGEList using filtered counts
DGE <- DGEList(counts = filtered, group = grouptype$group)

# Normalize data by read depth per sample
DGE <- calcNormFactors(DGE, method = "TMM")

# Obtain summary of DGEList
DGE

# Distance between samples based on gene expression
plotMDS(DGE, method = "bcv")
# Outliers? SUN_2_SILIQUE & IMB211_SHADE_3_LEAF --> not anymore when apical &floral meristem removed
plotMDS(DGE, method = "bcv", top=5)

# Boxplots - Added 1 to all filtered counts to avoid log2(0)
png("Final_normalized_boxplots.png")
boxplot(log2(cpm(DGE)+1),ylab="log2(Counts)",
        main="Normalized Filtered RNA-Seq Counts",xlab="Sample")
dev.off()
###################################
#Treatment only DE
# Make design matrix
design <- model.matrix(~trmt, data = grouptype)
rownames(design) <- colnames(DGE)

# Estimate dispersion
# a. Overall dispersion
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE)
# b. Trended dispersion
DGE <- estimateGLMTrendedDisp(DGE, design)
# c. Gene-wise dispersion
DGE <- estimateGLMTagwiseDisp(DGE, design)

# Obtain summary of DGEList
DGE

# Plot biological coefficient of variation
png("Trt_biological_coefficient_of_variation.png")
plotBCV(DGE)
dev.off()

# Fit model using design matrix
fit <- glmFit(DGE, design)

# Identify differentially expressed genes 
# between treatments using likelihood ratio test
lrt.trmt <- glmLRT(fit, coef = "trmtSHADE")

# Report topTags
topTags(lrt.trmt)

# Report how many are up- or down-regulated in shade relative to sun
DE <- decideTestsDGE(lrt.trmt, adjust.method = "fdr", p.value = 0.01)
summary(DE)

# Save results for all genes regardless of significance
results.trmt <- as.data.frame(topTags(lrt.trmt, n = length(rownames(DGE$counts))))
results.trmt<-cbind(results.trmt, V1=rownames(results.trmt))

# Write results to a tab-delimited text file
write.table(results.trmt,"Final_SHADE_SUN_results.txt",sep="\t",row.names=FALSE,quote=FALSE)

# Make MA plot
pdf("Final_MAplottrmt.pdf")
detags <- rownames(DGE)[as.logical(DE)]
plotSmear(lrt.trmt, de.tags = detags)
abline(h = c(-1, 1), col="blue")
dev.off()

#################################
#G & E
# Make design matrix
design <- model.matrix(~ geno + trmt, data = grouptype)
rownames(design) <- colnames(DGE)

# Estimate dispersion
# a. Overall dispersion
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE)
# b. Trended dispersion
DGE <- estimateGLMTrendedDisp(DGE, design)
# c. Gene-wise dispersion
DGE <- estimateGLMTagwiseDisp(DGE, design)

# Obtain summary of DGEList
DGE

# Plot biological coefficient of variation
png("Final_biological_coefficient_of_variation.png")
plotBCV(DGE)
dev.off()

# Fit model using design matrix
fit <- glmFit(DGE, design)

# Identify differentially expressed genes 
# between genotypes using likelihood ratio test
lrt.geno <- glmLRT(fit, coef = "genoR500")

# Report topTags
topTags(lrt.geno)

# Report how many are up- or down-regulated in R500 relative
# to IMB211
DE <- decideTestsDGE(lrt.geno, adjust.method = "fdr", p.value = 0.01)
summary(DE)

# Save results for all genes regardless of significance
results.geno <- as.data.frame(topTags(lrt.geno, n = length(rownames(DGE$counts))))

# Write results to a tab-delimited text file
write.table(results.geno,"Final_IMB211_R500_results.txt",sep="\t",row.names=FALSE,quote=FALSE)

# Make MA plot
pdf("Final_MAplotgeno.pdf")
detags <- rownames(DGE)[as.logical(DE)]
plotSmear(lrt.geno, de.tags = detags)
abline(h = c(-1, 1), col="blue")
dev.off()

########################
plotDE <- function(genes, dge, sample.description) {
  require(ggplot2)
  require(reshape2)
  tmp.data <- t(log2(cpm(dge[genes,])+1))
  tmp.data <- merge(tmp.data,grouptype,by.x="row.names",by.y="sample")
  tmp.data <- melt(tmp.data,value.name="log2_cpm",variable.name="gene")
  pl <- ggplot(tmp.data,aes(x=geno,y=log2_cpm,fill=trmt))
  pl <- pl + facet_wrap( ~ gene)
  pl <- pl + ylab("log2(cpm)") + xlab("genotype")
  pl <- pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}

plotDE(rownames(results.geno),DGE,grouptype)

plotDE("Bra009785",DGE,grouptype)
################################

