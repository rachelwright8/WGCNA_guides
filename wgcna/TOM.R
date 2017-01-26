library(WGCNA)
options(stringsAsFactors = FALSE)

lnames = load(file="SamplesAndTraits.RData")

softPower=20 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

save(TOM, dissTOM, file = "TOMdissTOM.RData")

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
pdf(file="dendrogram_thresh9_signed.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()

minModuleSize=30 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)

dynamicColors= labels2colors(dynamicMods)

MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 20)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_nomerge.RData")

sizeGrWindow(7,6)
pdf(file="nomerge.pdf", width=20, height=20)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()

MEDissThres= 0.0
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="nomerge.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_nomerge.RData")

