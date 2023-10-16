library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)


data<-read.csv("all_des_output.csv")
clean_data<-data[complete.cases(data),]
up_genes<-clean_data$genename[(clean_data$padj<=0.05)&(clean_data$log2FoldChange>1)]

ensLookup <- up_genes
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize = 10, maxGSSize = 500
                )
head(summary(ego))
pdf("RNA up GO.pdf")
dotplot(ego, showCategory=10,title="Enrichment GO(Up)")
dev.off()
write.csv(summary(ego),"enrich-Up.csv",row.names =FALSE)

down_genes<-clean_data$genename[(clean_data$padj<=0.05)&(clean_data$log2FoldChange< -1)]
ensLookup<-down_genes
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize = 10, maxGSSize = 500
                )
head(summary(ego))
pdf("RNA down GO.pdf")
dotplot(ego, showCategory=10,title="Enrichment GO(Down)")
dev.off()
write.csv(summary(ego),"enrich-down.csv",row.names =FALSE)

