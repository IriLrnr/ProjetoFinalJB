library("tidyr")

# Here a vector with the file names is created 
file.names<-dir("./data/position/")

# A new table is made, with all the times
ind.loc <- data.frame()

for (i in 1:length(file.names)){
    data.loc <- read.csv(paste("./data/position/", file.names[i], sep=""), sep=";", head=FALSE)
    ind.loc <- rbind(ind.loc, data.loc)
}

colnames(ind.loc)=c("id","x","y","sp","gen")

# Export the final table
write.csv(x = ind.loc, 
          file = "./data/indloc.csv", 
          row.names = FALSE)

