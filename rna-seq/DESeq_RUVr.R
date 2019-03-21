# RUV RNA-seq pipeline 0.1

library(RUVSeq)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(pamr)

# Retrieve gene quantification matrix
files <- list.files("female") #retrieve all files in female folder.
# feature Count
x <- readDGE(files, path="female", columns=c(1,2), sep=" ") #readDGE Reads and merges a set of text files containing gene expression counts.
class(x)
dim(x)

# Save R object as output (rda)
save(x, file="raw_rnaseq_subset.rda")
# save raw matrix into 

samplenames <- substring(colnames(x), 0, 23)
samplenames
colnames(x) <- samplenames
total_samples = 6

## ----filter lowly expressed genes ##
cpm <- cpm(x)
table(rowSums(x$counts>0)==total_samples)
keep.exprs <- rowSums(cpm>3)>=2
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----create sample group for DE analysis---- ##
expo <- as.factor(rep(c("Control","Arsenic"), c(3,3)))
x$samples$expo <- expo
group <- factor(paste(expo, sep="."))
lcpm <- cpm(x, log=FALSE)
genes <- rownames(x)


## ----store_data-------------------------------------------
set <- newSeqExpressionSet(as.matrix(x),
                           phenoData = data.frame(expo, row.names=colnames(x)))



## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
library(RColorBrewer)
par(mfrow=c(1,2))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[expo], cex=0.8)
plotPCA(set, col=colors[expo], cex=0.8)

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
design <- model.matrix(~expo, data=pData(set))
y <- DGEList(counts=counts(set), group=expo)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
# RUVr needs res values using deviance
res <- residuals(fit, type="deviance")

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05, lfc=0.583))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
#plotMD(lrt, pch=16, cex=1.2)
abline(h=c(-0.5, 0.5), col="blue")

## ----emp_RUVr, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
set4 <- RUVr(set, genes, k=3, res) #RUVr
par(mfrow=c(1,2))
pData(set4)
plotRLE(set4, outline=FALSE, ylim=c(-4, 4), col=colors[expo])
plotPCA(set4, col=colors[expo], cex=0.8)

## DESeq2 implemantation
cts = counts(set4)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = pData(set4),
                              design = ~ W_1 + W_2 + W_3 + expo)
dds <- DESeq(dds)
res <- results(dds)
res

sum(res$padj < 0.001, na.rm=TRUE)

design <- model.matrix(~expo + W_1, data=pData(set4))
y <- DGEList(counts=counts(set4), group=expo)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#  
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05, lfc=1.0))

par(mfrow=c(1,2))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
plotMD(lrt, pch=16, cex=0.5)
abline(h=c(-0.5, 0.5), col="blue")

tab <- topTags(lrt, n=2000, sort.by="p.value", p.value=0.05)
tab_all <- topTags(lrt, n=20000, sort.by="p.value", p.value=1)
write.table(tab_all, file="promoter_Male_RUVr_all.txt")
write.table(tab, file="promoter_Male_edgeR_RUVr.txt")
write.table(lcpm, file="promoterMale_RUVr_list_lcpm.txt")
########## edgeR ##########

########## data visualization ##############
library(gplots)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
genelist<-which(de[,1] == 1 | de[,1] == -1)
#genelist<-which(de[,1] == 1)
sub_set = lcpm[genelist, 1:total_samples]
col.cell <- c("blue","red")[expo]
col.cell = col.cell[1:total_samples]

# Plot the heatmap
heatmap.2(sub_set,col=rev(morecols(50)),trace="none", 
          main="Differentially expressed genes",
          ColSideColors=col.cell,scale="row",
          Rowv=TRUE,Colv=TRUE)
############################################

