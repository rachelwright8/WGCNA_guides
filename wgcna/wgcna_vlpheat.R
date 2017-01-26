setwd("~/Documents/virusrna/DESeq2LRT/wgcna/")
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

# #Load expression data-----
# dis0 = read.csv("../4August2016_aims2013_vlp_heat_BM3_RLDandPVALS.csv")
# names(dis0)
# dis = dis0[c(2:20)]
# row.names(dis) = dis0$X
# dis$X = NULL
# datExpr0 = as.data.frame(t(dis))
# row.names(datExpr0) = sub("X", "", row.names(datExpr0))
# 
# 
# gsg = goodSamplesGenes(datExpr0, verbose = 3)
# gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts
# 
# #Load trait data-----
# d = row.names(datExpr0)
# dam1 = substr(d,1,1)
# sire1 = substr(d,2,2)
# vlp1 = substr(d,3,3)
# heat1 = substr(d,4,4)
# 
# traits = as.data.frame(cbind(d, dam1, sire1, vlp1, heat1))
# traits
# 
# traits$dam = ifelse(traits$dam1=="1","1","0")
# traits$sire = ifelse(traits$sire1=="A","1", ifelse(traits$sire1=="C","2","3"))
# traits$vlp = ifelse(traits$vlp1=="V","1","0")
# traits$heat = ifelse(traits$heat1=="H","1","0")
# 
# datTraits0 = traits[c(1,6:9)]
# datTraits0
# 
# # read in trait data
# vegaData = read.csv("traits4WGCNA.csv")
# head(vegaData)
# 
# datTraits = merge(datTraits0, vegaData, by=1)
# head(datTraits)
# row.names(datTraits) = datTraits$d
# datTraits$d = NULL
# 
# # check that samples are lined up properly in Expression and Trait data
# match(row.names(datExpr0), row.names(datTraits)) == c(1:19)
# 
# save(datExpr0, datTraits ,file = "SamplesAndTraits.RData")

# Load Data---------------------

load("SamplesAndTraits.RData")
load("Network_nomerge.RData")
load("Network_signed_nomerge.RData")


# START HERE AFTER TACC -------------------------

# Make Modules ---------------------------
softPower = 14
dynamicColors= labels2colors(dynamicMods)

quartz()
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

# Merge modules whose expression profiles are very similar -------------

MEList = moduleEigengenes(datExpr0, colors= dynamicColors)
MEs = MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

quartz()
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.4
abline(h=MEDissThres, col="red")

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors
mergedMEs = merge$newMEs

quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= paste("SamplesAndColors_soft",softPower,"merge",MEDissThres,"_signed.RData",sep=""))
                                                            
# load ME data ------------------------------------                                                            
#load("SamplesAndColors_soft9merge0.5.RData")

# organize traits
head(datTraits)
datT = datTraits
datT$vlp.heat = ifelse(substr(row.names(datT),3,4)=="VH",1,0)
head(datT)
names(datT)
datT = datT[c(1:4,8,5:7)]

# Relate modules to traits---------------------
datt = datExpr0

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#Correlations of genes with eigengenes
moduleGeneCor = cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);
moduleTraitCor = cor(MEs, datT, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Module-trait heatmap ----------
quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datT),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Gene significance by Module membership scatterplots-------------------
whichTrait="vlp" #Replace this with the trait of interest

quartz()
nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(datTraits[,whichTrait]);
names(selTrait) = whichTrait
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
par(mfrow=c(3,4))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>12) {
    quartz()
    par(mfrow=c(3,4))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = module,mgp=c(2.3,1,0))
}

# Eigengene heatmap--------
which.module="antiquewhite4" #replace with module of interest
datME=MEs
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr0[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=c(row.names(datt)), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample", las=2)

#Heatmap for top-kME in a module with gene names -----------------
rld = read.csv("../4August2016_aims2013_vlp_heat_BM3_RLDandPVALS.csv")
names(rld)

a.rld = rld[c(2:20)] #Columns with rlog data
row.names(a.rld) = rld$X
head(a.rld)

allkME = as.data.frame(signedKME(t(a.rld), MEs))

gg = read.delim("~/Documents/genomes/amil_moya/amil_iso2gene.tab", sep="\t", header=F)
head(gg)
library(pheatmap)

whichModule = "antiquewhite4"
top = 50

modcol = paste("kME",whichModule,sep="")
sorted = a.rld[order(allkME[,modcol],decreasing=T),]
hubs = sorted[1:top,]
# attaching gene names
gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

names(hubs) = sub("X", "", names(hubs))

quartz()
pheatmap(hubs,scale = "row",col=contrasting,border_color = NA, main = paste(whichModule,"top",top,"kME",sep = ""), cluster_cols = T)

#Fisher of Module vs Whole Dataset for GO Analysis -------------
datME = moduleEigengenes(datExpr0,moduleColors)$eigengenes
datKME = signedKME(datExpr0, datME, outputColumnName="MM.")
genes = names(datExpr0)
geneInfo0 = data.frame(gene=genes,moduleColor=moduleColors)
color = data.frame(geneInfo0,datKME) #these are from your original WGCNA analysis 
head(color)

test = allkME
names(test)
#change test column to first empty column number
for (i in row.names(test)) { 
  gi=color[color$gene==i,]
  if (length(gi[,1])>0){
    test[i,25]<-gi$moduleColor    
  }    
}

test$V25 = as.factor(test$V25)
test[,25] = as.factor(test[,25])

for(i in names(test[c(1:24)])){
  cat <- test
  cat$fisher <- ifelse(cat$V25==gsub("kME","",i),1,0) #make column called "fisher"
  catt <- cat[c(26)] # subset for column called "fisher"
  write.csv(catt,file=paste("fisher_",gsub("kME","",i),".csv",sep=""),quote=F,row.names=T)
}


#RLD by module for heatmaps, PCA, etc ------------------
moduleColorList <-names(MEs)
moduleColorList <- sub("ME","",moduleColorList)
for (i in moduleColorList) { 
  cands <- names(datExpr0[moduleColors==i])
  c.rld <- rld[rld$X %in% cands,]
  write.csv(c.rld,paste("rld_",i,".csv",sep=""),quote=F, row.names=F)
}

# GO by kME within module only ---------
allkME  <-  as.data.frame(signedKME(t(a.rld), MEs))
head(allkME)

for (i in moduleColorList) { 
  modColName=paste("kME",i,sep="")
  modkME=as.data.frame(allkME[,modColName])
  modkME[,"gene"]=row.names(allkME)
  names(modkME)=c(modColName,"gene")
  modkME=modkME[c(2,1)]
#--subset for only genes within the module
  cands=names(datExpr0[moduleColors==i])
  inmodkME=modkME[(modkME$gene %in% cands),]
  write.csv(inmodkME,file=paste("kME",i,".csv",sep=""),quote=F, row.names=F)
}

# GO by kME + fisher (0 for not-in-module) -----------

for (i in moduleColorList) { 
  modColName=paste("kME",i,sep="")
  modkME=as.data.frame(allkME[,modColName])
  modkME[,"gene"]=row.names(allkME)
  names(modkME)=c(modColName,"gene")
  modkME=modkME[c(2,1)]
  #--subset for only genes within the module
  cands=names(datExpr0[moduleColors==i])
  inmodkME=modkME[(modkME$gene %in% cands),]
  outmodkME=modkME[!(modkME$gene %in% cands),]
  outmodkME[2]=c(rep(0,length(outmodkME[,2])))
  inmodkME=rbind(inmodkME,outmodkME)
  write.csv(inmodkME,file=paste("kME_fisher",i,".csv",sep=""),quote=F, row.names=F)
}

