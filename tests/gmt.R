library(BioQC)

gsItem1 <- gsitem("GeneSet1", "Desc", c("G1", "G2", "G6"))
gsItem2 <- gsitem("GeneSet2", "Desc 2", c("G2", "G3", "G4", "G5", "G6", "G7"))
gsItem3 <- gsitem("GeneSet3", "Desc 3", c("G4", "G5", "G6"))
gsItem4 <- gsitem("GeneSet4", "Desc 4", "G6")

gss <- list(gsItem1, gsItem2, gsItem3, gsItem4)
gsList <- as(gss, "gslist")

show(gsItem1)
show(gsItem2)
show(gsItem3)
show(gsItem4)
show(gsList)

gsGenes(gsItem1)
gsGenes(gsItem2)
gsGenes(gsItem3)
gsGenes(gsItem4)
gsGenes(gsList)
gsUniqGenes(gsList)

(gsDf <- as.data.frame(gsList))
as(gsList, "data.frame")

(gsMat <- as.matrix(gsList))
as(gsMat, "matrix")

tfidf(gsList)
gsListTFIDF <- gsTfIdf(gsList)
gsWeights(gsListTFIDF)

gsList2 <- gsList
print(gsList2)
gsWeights(gsList2) <- list(2,2,2,2)
gsWeights(gsList2)
gsWeights(gsList3 <- gsClearWeights(gsList2))

gsHasWeights(gsList)
gsHasWeights(gsList2)
gsHasWeights(gsList3)
gsHasWeights(gsListTFIDF)
