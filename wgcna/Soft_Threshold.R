source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
biocLite("WGCNA")

library(WGCNA)
options(stringsAsFactors = FALSE)

lnames= load(file="SamplesAndTraits.RData")

dim(datExpr0)
dim(datTraits)
powers= c(c(1:10), seq(from =12, to=40, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
pdf(file="soft_threshold_signed.pdf", width=30, height=15)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.80, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
