library(DESeq2)
sampleNames <- c("ctrl_1", "ctrl_2","ctrl_3", "iko_1", "iko_2","iko_3")
data <- read.table("out.bed", header=TRUE, quote="\t", skip=1)
# 前六列分别是Geneid	Chr	Start	End	Strand	Length
# 我们要的是count数，所以从第七列开始
names(data)[7:12] <- sampleNames
countData <- as.matrix(data[7:12])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("ctrl", "ctrl","ctrl", "iko", "iko","iko"),type=c("1","2","2","1","2","2"))
rownames(database) <- sampleNames
database$type <- as.factor(database$type)

## 设置分组信息并构建dds对象
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition+type)
dds <- dds[ rowSums(counts(dds)) > 1, ]

## 使用DESeq函数估计离散度，然后差异分析获得res对象
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "res_des_output.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
promoter<-read.table("/lustre/home/fengyang/Ref/mm10/mm10_promoter_1k_1k.bed")
resdata$genename="None"
for(i in 1:length(resdata$Row.names)){
	temp_gene=promoter$V5[promoter$V4==resdata$Row.names[i]]
	if(length(temp_gene)>0){
		resdata$genename[i]=temp_gene
	}
}
write.csv(resdata, "all_des_output.csv", row.names=FALSE)

# library(DESeq2)
pdf("plot.pdf")
plotMA(res, main="DESeq2", ylim=c(-2, 2))
dev.off()


library(ggplot2)

# 这里的resdata也可以用res_des_output.csv这个结果重新导入哦。
# 现在就是用的前面做DESeq的时候的resdata。
resdata$change <- as.factor(
	ifelse(
		resdata$padj<0.05 & abs(resdata$log2FoldChange)>1,
		ifelse(resdata$log2FoldChange>1, "Up", "Down"),
		"NoDiff"
	)
)
resdata<-resdata[complete.cases(resdata),]

valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
	geom_point(alpha=0.8, size=1) + 
	theme_bw(base_size=10) + 
	theme(
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank()
	) + 
	ggtitle("DESeq2 Valcano") + 
	# scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
	geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
	geom_hline(yintercept=-log10(0.05), lty=1, col="gray", lwd=0.5)+xlim(-10,10)
pdf("valcano.pdf")
valcano
dev.off()


# library(ggplot2)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape = type)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text(aes(label=name),vjust=2)
pdf("PCA.pdf")
pca
dev.off()

# rld <- rlog(dds)
# pcaData <- plotPCA(rld, intgroup=c("condition", "name"), returnData=T)
# percentVar <- round(100*attr(pcaData, "percentVar"))
# pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=name)) + 
# 	geom_point(size=3) + 
# 	ggtitle("DESeq2 PCA") + 
# 	xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
# 	ylab(paste0("PC2: ", percentVar[2], "% variance"))


assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$type)

pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape = type)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text(aes(label=name),vjust=2)
pdf("PCA_BatchEffectRemoved.pdf")
pca
dev.off()
write.table(assay(vsd),file='pasilla_batchCorrectedVSD.txt',sep="\t",row.names=TRUE,col.names=TRUE)
# library(pheatmap)

# # select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
# # nt <- normTransform(dds)
# # log2.norm.counts <- assay(nt)[select,]
# # df <- as.data.frame(colData(dds)[, c("name", "condition")])
# # pdf("pheatmap.pdf")
# # pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
# # dev.off()

# promoter<-read.table("/lustre/fengyang/Xiaona/SC_G4/G4/G4_p001_peaks/ASC48_vs96unique_promoter.bed",sep="\t")
# # promoter<-read.table("/lustre/fengyang/Xiaona/SC_G4/G4/G4_p001_peaks/ASC96_vs48unique_promoter.bed",sep="\t")
# nt <- normTransform(dds)
# select <- rownames(nt)[rownames(nt) %in% promoter$V4]
# log2.norm.counts <- assay(nt)[select,]
# df <- as.data.frame(colData(dds)[, c("name", "condition")])
# pdf("G4_pheatmap.pdf")
# pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
# dev.off()


# temp<-resdata[resdata$Row.names %in% promoter$V4,]
# temp<-temp[complete.cases(temp),]

# temp_up<-temp[temp$change=="Up",]
# specific<-read.table("muscle_pescific.bed",sep="\t")
# intersect(specific$V1,temp_up$genename)