library("tidyr")

setwd("~/ProjetoFinalJB/data/position") #Bad

# Here a vector with the file names is created 
file.names <- dir()

# A new table is made, with all the times
ind.loc <- data.frame()

for (i in 1:length(file.names)){
  data.loc <- read.csv(paste(file.names[i]), head=FALSE, sep=";")
  ind.loc <- rbind(ind.loc, data.loc)
}

colnames(ind.loc)=c("id","x","y","sp","time")

# Export the final table
write.csv(x = ind.loc, 
          file = "~/ProjetoFinalJB/data/indloc.csv", 
          row.names = FALSE)


