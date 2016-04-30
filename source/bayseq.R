source("https://bioconductor.org/biocLite.R")
biocLite("baySeq")
library("baySeq", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")

# Read in data file
dat <- read.table(file = "/Users/mblai/Documents/BigData/Brassica_parents/GC_merged_v1.5_mapping.tsv",
                  header = T)

# Remove * row and make gene IDs row labels
data <- dat[2:41921, 2:67]
rownames(data) <- dat[2:41921,1]
# Shorten the names
colnames(data)<-sub(".1.merged.bam", "", colnames(data))
#remove apical meristem
data<-data[,c(-1,-7,-32,-39,-51)]
#remove floral meristem
data<-data[,c(-30,-36)]
# Convert NAs to zeros
data[is.na(data)]<-0
# Round decimals up to nearest integer
data<-ceiling(data)

# Filter out lowly expressed genes
# filter genes > 10 reads in >= 3 samples
filtered <- data[rowSums(data > 10) >= 3,]
genes<-as.matrix(filtered)

#set replicate structure for genotype by environment
#expect some equivalent expression between replicates
replicates<-c(rep("IMB211_SHADE",14),rep("IMB211_SUN",15),
              rep("R500_SHADE",15),rep("R500_SUN",15))

#identify predictions of no differential expression & differential
  #expression
groups<-list(NDE=c(rep(1,14),rep(2,15),rep(3,15),rep(4,15)),
      DE= c(1,2,3,4,5,1,2,3,4,5,2,3,4,5,rep(6:10,3),
          rep(11:15,3), rep(16:20,3)))

#combine count data & groups in countData object
CD <- new("countData", data = genes, 
          replicates = replicates, groups = groups)

#infer library sizes from the data
libsizes(CD)<-getLibsizes(CD)

CD@annotation <- data.frame(rownames(filtered))

CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
CD@estProps
CD@posteriors[1:10,]










#MA plot
plotMA.CD(CD, samplesA = c("IMB211_SHADE","IMB211_SUN"),
          samplesB = c("R500_SHADE","R500_SUN" ))

plotMA.CD(CD, samplesA ="IMB211_SHADE" , samplesB ="R500_SHADE" )
######################################
#SHADE VS SUN ONLY
#######################################################
#replicates<-c(rep("SHADE",14),rep("SUN",15),
 #             rep("SHADE",15),rep("SUN",15))

#identify predictions of no differential expression & differential
#expression
#groups<-list(NDE=rep(1,59), DE=c(rep(1,14),rep(2,15),rep(1,15),rep(2,15)))
#groups<-list(NDE=c(rep(1,14),rep(2,15),rep(3,15),rep(4,15)),
#             DE= c(1,2,3,4,5,1,2,3,4,5,2,3,4,5,rep(6:10,3),
 #                  rep(11:15,3), rep(16:20,3)), 
  #           annotation=c(rep("IMB211_SHADE",14),
   #                       rep("IMB211_SUN",15),
    #                      rep("R500_SHADE",15),rep("R500_SUN",15))
#)

#combine count data & groups in countData object
#CD <- new("countData", data = genes, 
 #         replicates = replicates, groups = groups)

#infer library sizes from the data
#libsizes(CD)<-getLibsizes(CD)

#MA plot
#plotMA.CD(CD, samplesA = c("IMB211_SHADE","IMB211_SUN"),
 #         samplesB = c("R500_SHADE","R500_SUN" ))

#plotMA.CD(CD, samplesA ="SHADE" , samplesB ="SUN", 
#          )

















