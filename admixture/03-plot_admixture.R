#largely based off of tutorial at https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html

#load input files
tbl <- read.table("67_auto.2.Q") #Q file from ADMIXTURE program
indTable <- read.csv("indTable.csv", header=T) #sample and pop labels

#plotting raw Q file without pop label adjustments
barplot(t(as.matrix(tbl)), col=c("tomato2", "skyblue3"), xlab="Individual #", ylab="Ancestry", border=NA)

#working with pop labels to order the admixture map
mergedAdmixTable <- cbind(tbl, indTable)
ordered <- mergedAdmixTable[order(mergedAdmixTable$Pop),]

#instead of alphabetical order, changing the order of pops to be by latitude (matching fst plots, etc). Subsetting data to merge in custom order
AB <- subset(ordered, Pop=="AB")
FLD <- subset(ordered, Pop=="FLD")
TX <- subset(ordered, Pop=="TX")
PA <- subset(ordered, Pop=="PA")
re_ordered <- rbind(AB, PA, TX, FLD)

#label by individual
barplot(t(as.matrix(re_ordered)), col=c("tomato2", "skyblue3"), xlab="Individual #", ylab="Ancestry", border=NA, main="admixture mapping K=2")

#barNaming function for labeling plot by pop
barNaming <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#label by population
barplot(t(as.matrix(re_ordered)), col=c("tomato2", "skyblue3"), xlab="Population", ylab="Ancestry", border=NA, names.arg=barNaming(re_ordered$Pop), main="admixture mapping K=2")
