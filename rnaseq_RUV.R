# RUV

library(RUVSeq)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(pamr)

files <- list.files("./Multi-omics/featureCounts_pairwise")

x <- readDGE(files, path="./Multi-omics/featureCounts_pairwise", columns=c(1,3))
class(x)
dim(x)

save(x, file="./Multi-omics/raw_rnaseq_subset.rda")
# save raw matrix into 

samplenames <- substring(colnames(x), 0, 5)
samplenames
colnames(x) <- samplenames

gender <- as.factor(rep(c("Male"), c(8)))
expo <- as.factor(rep(c("FA","PM"), c(4,4)))

x$samples$gender <- gender
x$samples$expo <- expo
group <- factor(paste(gender, expo, sep="."))

## ----filter lowly expressed genes ##
cpm <- cpm(x)
table(rowSums(x$counts==0)==8)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
lcpm <- cpm(x, log=TRUE)
## ----store_data-------------------------------------------
expo <- as.factor(rep(c("FA","PM"), c(4,4)))
set <- newSeqExpressionSet(as.matrix(x),
                           phenoData = data.frame(expo, row.names=colnames(x)))
set

## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
library(RColorBrewer)
par(mfrow=c(1,2))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[expo])
plotPCA(set, col=colors[expo], cex=1.2)

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
design <- model.matrix(~expo, data=pData(set))
y <- DGEList(counts=counts(set), group=expo)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

## ----emp_ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set2 <- RUVg(set, empirical, k=1)
par(mfrow=c(1,2))
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[expo])
plotPCA(set2, col=colors[expo], cex=1.2)

design <- model.matrix(~expo + W_1, data=pData(set2))
y <- DGEList(counts=counts(set2), group=expo)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#  
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.01, lfc=0.5))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
#plotMD(lrt, pch=16, cex=1.2)
abline(h=c(-0.5, 0.5), col="blue")

tab <- topTags(lrt, n=2000, sort.by="p.value", p.value=0.01)
write.table(tab, file="mygenelist_edgeR.txt")
########## edgeR ##########

########## data visualization ##############
library(gplots)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
#genelist<-which(de[,1] == 1 | de[,1] == -1)
genelist<-which(de[,1] == 1)
sub_set2 = lcpm[genelist, 1:8]
col.cell <- c("blue","red")[sub_set2$samples$expo]
col.cell = col.cell[1:8]

# Plot the heatmap
heatmap.2(sub_set2,col=rev(morecols(50)),trace="none", 
          main="Differentially expressed genes",
          ColSideColors=col.cell,scale="row",
          Rowv=TRUE,Colv=TRUE)
############################################

