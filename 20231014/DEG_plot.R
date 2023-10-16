WT_RNA<-read.table("/lustre/fengyang/JunB/SC_RNA/ctrl/02quantification/genes.fpkm_tracking",header=T)
KO_RNA<-read.table("/lustre/fengyang/JunB/SC_RNA/iko/02quantification/genes.fpkm_tracking",header=T)
# 470246.1 WT
# 604257.3 KO
KO_RNA$FPKM<-KO_RNA$FPKM*(604257.3/470246.1)


for(i in 1:55240){
	index <- which((KO_RNA$gene_id==WT_RNA$gene_id[i])&(KO_RNA$locus==WT_RNA$locus[i]))

	if(length(index) == 1){
		WT_RNA$FPKM_KO[i] <- KO_RNA$FPKM[index]
		WT_RNA$fold_change[i] <- (KO_RNA$FPKM[index]+1)/(WT_RNA$FPKM[i]+1)
		WT_RNA$fold_change_log[i] <- log2(WT_RNA$fold_change[i])
	}
	else{
		print(KO_RNA$gene_short_name[i])
	}
}

cut_off_logFC = 1           #差异倍数值

WT_RNA$change = ifelse(abs(WT_RNA$fold_change_log) > cut_off_logFC, 
                          ifelse(WT_RNA$fold_change_log> cut_off_logFC ,'Up','Down'),
                          'Stable')
result <- WT_RNA
sum(result$change=='Down') # 1514

sum(result$change=='Up') # 2669

library(ggplot2)
result_filter <- result[((result$FPKM_KO>0.5)|(result$FPKM)>0.5),]
#16400
p <- ggplot(data = result_filter, 
			mapping = aes(x = log2(FPKM+1), y = log2(FPKM_KO+1),color = change))+ ylim(0,12) + xlim(0,12) +
		theme(legend.title = element_text(face='bold'),axis.title.y = element_text(face='bold',size=15),axis.title.x = element_text(face='bold',size=15),plot.title = element_text(face='bold',size=20,hjust = 0.5)) + 
		ggtitle("Gene expression") +
  		xlab("log2(FPKM_WT+1)") + ylab("log2(FPKM_KO+1)")+
  		geom_point(alpha=0.8,size = 0.3)+
  		# geom_abline(intercept=1,slope=1)+geom_abline(intercept=-1,slope=1)
pdf("Gene_express.pdf")
p
dev.off()
sum(result_filter$change=='Down') # 1514

sum(result_filter$change=='Up') # 2669

write.csv(result_filter,file='result_filter.csv')


library(org.Mm.eg.db)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)


up <- result_filter[result_filter$change=="Up",]
down <- result_filter[result_filter$change=="Down",]
ensLookup<-up$gene_short_name
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
pdf("EnrichmentGO_UP.pdf")
dotplot(ego, font.size=10,title="Enrichment GO(UP)")
# barplot(ego, showCategory=20,title="EnrichmentGO_UP")
dev.off()
write.csv(summary(ego),"enrich-Up.csv",row.names =FALSE)


ensLookup<-down$gene_short_name
ego <- enrichGO(gene          = ensLookup,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
pdf("EnrichmentGO_down.pdf")
dotplot(ego, font.size=10,title="EnrichmentGO_down")
dev.off()
write.csv(summary(ego),"enrich-down.csv",row.names =FALSE)


over_gene_up<-intersect(up$gene_short_name,genes)
over_gene_down<-intersect(down$gene_short_name,genes)
pdf("over_gene_up_GO.pdf")
dotplot(ego, font.size=10,title="EnrichmentGO_down")
dev.off()
write.csv(summary(ego),"over_gene_up_GO.csv",row.names =FALSE)