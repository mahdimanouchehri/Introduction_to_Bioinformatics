setwd("H:/sharif/intro_bio/project/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
##source("GEOpatch.R")
series <- "GSE48558"
platform <- "GPL6244"

gset <- getGEO("GSE48558", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]





gr <- c(rep("AML",13),"Granulocytes","Granulocytes","B_Cells","T_Cells","Granulocytes","Granulocytes",
        rep("Monocytes",2),"B_Cells","T_Cells",rep("T_Cells",2),
        rep("T_Cells",2), "B_Cells","T_Cells","B_Cells","T_Cells",
        "CD34","CD34","CD34",rep("Granulocytes",7),rep("AML",2)
        ,"T_Cells",rep("AML",3),rep("B_Cells",7),"T_Cells",rep("Monocytes",4),"Granulocytes",rep("T_Cells",7))




fvarLabels(gset) <- make.names(fvarLabels(gset))
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
               "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
               "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
               "11111111111111111111")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
ex <- exprs(gset)


pdf("Results/boxplot3.pdf",width = 64,)
boxplot(ex)
dev.off()




####normalize

#ex <- normalizeQuantiles(ex)
#ex <- exprs(gset)

#pdf("Results/boxplotnormal.pdf",width = 64,)
#boxplot(ex)
#dev.off()





#pca
pc <- prcomp(ex)
pdf("Results/pc5.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()


ex.scale <- t(scale(t(ex),scale = F))
pc <- prcomp(ex.scale)
pdf("Results/pc_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()




pcr <- data.frame(pc$r[,1:3],Group=gr)
pdf("Results/pca_samples1.pdf")
ggplot(pcr,aes(PC1, PC2,color=Group)) + geom_point(size=3) + theme_bw()
dev.off()





#heatmap
pdf("Results/CorHeatmap3.pdf",width = 20,height = 20)
pheatmap(cor(ex),labels_row = gr,labels_col = gr)
dev.off()



####differential expression analysis
gr <- factor(gr)
gset$group <- gr
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gr)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
#cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(AML-CD34, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","P.Value","logFC"))

write.table(tT,"Results/AML__CD34.txt", row.names=F, sep="\t",quote=F)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
aml.up.genes <- sub("///.*","",aml.up.genes)
write.table(aml.up.genes, file = "Results/AML_CD34_UP.txt",quote=F,row.names = F, col.names = F)

aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
aml.down.genes <- sub("///.*","",aml.down.genes)
write.table(aml.down.genes, file = "Results/AML_CD34_down.txt",quote=F,row.names = F, col.names = F)