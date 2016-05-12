if(require("parallel")) cl <- makeCluster(8) else cl <- NULL


data(simData)
replicates <- c("simA", "simA", "simA", "simA", "simA",
                     "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1),
             DE = c(1,1,1,1,1,2,2,2,2,2))

CD <- new("countData", data = simData, replicates = replicates, groups = groups)
libsizes(CD) <- getLibsizes(CD)
plotMA.CD(CD, samplesA = "simA", samplesB = "simB",
     col = c(rep("red", 100), rep("black", 900)))
CD@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))

CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
CD@estProps
CD@posteriors[1:10,]

topCounts(CD, group = "DE")
