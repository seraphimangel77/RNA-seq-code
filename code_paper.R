# heatmap -----------------------------------------------------------------

sta.data <- read.table("/home/data/t050204/Data/ZYS/A_B.txt", header=T,sep="\t")

library(pheatmap)
DE_A1=c("AMOTL2","ANKRD1","ATAD2","CEP57","CTGF","CYR61","GADD45B","MYBL1","MYCMYC")
YAP=sta.data[DE_A1,]
YAP=YAP[,-2]
DE_A2=c("PANK2","POLA2","PSMC3IP","RAD18","RBM22","SAMD4A","TIMELESS","USP36","ZNHIT6")
YAP=YAP[DE_A2,]
YAP=YAP[,-1]
pheatmap(YAP, 
         #annotation_row=info, 
         #annotation_col=info, 
         show_colnames = TRUE, 
         show_rownames=TRUE, 
         fontsize=5, # 字体大小
         color = colorRampPalette(c('#4450A1','#ffffff','#BC0a00'))(50),
         annotation_legend=TRUE, 
         border_color= NA, 
         scale="row",  
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         angle_col = 90,
         annotation_colors = ann_colors,
         cellwidth=10,
         cellheight=10,
         filename="YAP_heatmap_all.pdf"
)


# GSEA --------------------------------------------------------------------

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
raw=read.csv("count.csv",sep = ",",header=TRUE)
df=raw[,2:8]
df=df[!duplicated(df$gene_name),]
rownames(df)=df$gene_name
df=df[,-7]
rt=df
rt=as.matrix(rt)
exp = rt[,1:ncol(rt)] 
barcode=colnames(exp)
normal=df[,1:3]
tumor=df[,4:6]
barcode_normal=colnames(normal)
barcode_tumor=colnames(tumor)
dataSmTP <- barcode_tumor
dataSmNT <- barcode_normal
exp = exp[,c(dataSmNT,dataSmTP)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=c(rep("normal",3),rep("tumor",3))
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
keep <-rowSums(y$counts) > 50
y <- y[keep,,keep.lib.size=FALSE]
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor")) 
ordered_tags <- topTags(et, n=100000)


allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
newData=y$counts
write.csv(allDiff,file="case_vs_control_RNA_edgerOut.csv",quote=F)
write.csv(newData,file="RNA_normalizeExp_TCGA.csv",quote=F) 
etSig <- allDiff[which(allDiff$PValue < 0.05 & abs(allDiff$logFC) > 0),] 
etSig[which(etSig$logFC > 0), "up_down"] <- "Up"
etSig[which(etSig$logFC < 0), "up_down"] <- "Down"
write.csv(etSig,file="case_vs_control_P0.05_FC2_RNA_edgerOut.csv",quote = F)
write.csv(newData[rownames(etSig),],file="RNA_diffExp.csv",quote=F) 
allDiff$logP <- -log10(allDiff$FDR)
allDiff <- cbind(symbol=rownames(allDiff),allDiff)
allDiff$Group = "not-signigicant"
allDiff$Group[which((allDiff$FDR < 0.05) & (allDiff$logFC > 0.2))]="up-regulated" 
allDiff$Group[which((allDiff$FDR < 0.05) & (allDiff$logFC < -0.2))]="down-regulated"
allDiff$Label=""
allDiff <- allDiff[order(allDiff$FDR),]
upgenes <- head(allDiff$symbol[which(allDiff$Group=='up-regulated')],50000)
downgenes <- head(allDiff$symbol[which(allDiff$Group=='down-regulated')],50000)
DEGs<-c(upgenes,downgenes)
organism = 'hsa' 
OrgDb = 'org.Hs.eg.db'

df=DEGs
geneEntrezID <- bitr(df, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id=geneEntrezID$ENTREZID
kegg = enrichKEGG(id, 
                  organism = 'human', 
                  keyType = 'kegg',
                  #keyType='ENTREZID',
                  pvalueCutoff = 1,
                  pAdjustMethod = 'BH',
                  minGSSize = 10,
                  maxGSSize = 500,
                  qvalueCutoff = 1,
                  use_internal_data = FALSE)
dotplot(kegg, showCategory=10) +
  theme(axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.title = element_text(size = 20))
ggsave("KEGG_FC1_A1B1.pdf",width = 9,height = 12)


UP=etSig[which(etSig$up_down=="Up"),]
DOWN=etSig[which(etSig$up_down=="Down"),]

UP=UP[,-2]
UP=UP[,-3:-4]
DOWN=DOWN[,-2]
DOWN=DOWN[,-3:-4]

need_DEG=rbind(UP,DOWN)
need_DEG$SYMBOL <- rownames(need_DEG)
df <- bitr(rownames(need_DEG), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb)
need_DEG <- merge(need_DEG, df, by='SYMBOL')

df_all=need_DEG
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]

gene_fc = df_all_sort$logFC
names(gene_fc) <- df_all_sort$ENTREZID

KEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1)

KEGG_kk <- DOSE::setReadable(KEGG, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')

kk_gse <- KEGG_kk

gseaplot2(kk_gse,
          title = "Hippo signaling pathway",  
          "hsa04390", 
          color="#E35858", 
          base_size = 0, 
          subplots = 1:2, 
          pvalue_table = F) 

dev.off()
dev.new()




# Cox ---------------------------------------------------------------------

library(grid)
library(forestploter)
library(foreign)
library(survival)
library("survminer")

bc<-read.csv("group_cox.tsv",sep='\t',header=TRUE)
bc$HGS<-as.factor(bc$Group1)
bc$Time<-as.numeric(bc$Time)
bc$Survival<-as.numeric(bc$Survival)
bc$Age<-as.numeric(bc$Age)
bc$Sex<-as.factor(bc$Sex)
f1<- coxph(Surv(Time, Survival) ~ HGS+Age+Sex+Stage, data=bc)
ggforest(f1,data = bc)
ggsave("forest_G1.pdf",width = 7,height = 5,dpi=500)




