library(data.table)
bay<-fread(file.choose())
edge<-fread(file.choose())

setkey(bay, V1)
setkey(edge, V1)

sharedDE <- merge(edge, bay)
write.csv(sharedDE, "shared_DEgenes.csv")
