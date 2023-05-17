library(pheatmap)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(grDevices)
library(RColorBrewer)

themo_demo= theme_bw()+
  theme(
    text = element_text(size=10),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text=element_text(color='black'),
    plot.title = element_text(hjust = 0.5),
    title=element_text(size = 12),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10))

###################################### load all required data before RNA-seq data analysis ####################################
load("/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/CD4_TF_TG_analysis_Step1_2023.Rdata")

###############################################################################################################################
###############                                      RNA-seq data (DEG)                                           #############
###############################################################################################################################
load("/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/Fig4_RNA_DEseq2.by.group_Ab_PI_2023.RData") ## produced by "Fig4_RNA_DEseq2.by.group.Ab_PI.R"

trd <- log2(2)
UP_DEG.TCR <- subset(DEG.all_act.up_dataframe, S.240_vs_NC.0 > trd) 
UP_DEG.PI <- subset(DEG.all_act.up_dataframe, PIC.240_vs_NC.0> trd) 

y1 <- calculate.overlap(x= list(TCR = unique(UP_DEG.TCR$Symbol),
                                PI = unique(UP_DEG.PI$Symbol)))

coDEG.TCR <- subset(UP_DEG.TCR, UP_DEG.TCR$Symbol %in% coTGs$target)
coDEG.PI <- subset(UP_DEG.PI, UP_DEG.PI$Symbol %in% coTGs$target)

y2 <- calculate.overlap(x= list(TCR = unique(coDEG.TCR$Symbol),
                                PI = unique(coDEG.PI$Symbol)))

UP_DEGs <- subset(DEG.all_act.up_dataframe, PIC.240_vs_NC.0> trd | S.240_vs_NC.0> trd) ##DEG.all.up_dataframe
# UP_DEGs <- UP_DEGs[order(-UP_DEGs$S.240_vs_NC.0),]
UP_DEG.TCR.PI <- subset(UP_DEGs, UP_DEGs$Symbol %in% y1[["a3"]]) #  705row, 669 gene
coTGs.DEG.TCR.PI <- subset(UP_DEG.TCR.PI, UP_DEG.TCR.PI$Symbol %in% coTGs$target) #104 row, 96 gene

setwd("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/")
# write.csv(coTGs.DEG.TCR.PI, "./combinatory_effect/coTGs.DEG.TCR.PI.csv", row.names = F)


######## from "/storage/labnas5/HuangWen2023/T cell dynamics paper/"
AUC.TCR.PI.2023 <- read.csv("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/For_Fig4.AUC_mean_TCR_PI_2023.csv")
AUC.CAR.2023 <- read.csv("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/For.Fig4.AUC_mean_CAR_2023.csv")
AUC.2023 <- rbind(AUC.TCR.PI.2023, AUC.CAR.2023)

# AUC.TCR.PI.2023 <- read.csv("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/Fig4.AUC_mean_TCR_PI.csv")
# AUC.CAR.2023 <- read.csv("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/Fig4.AUC_mean_CAR.csv")
# AUC.2023 <- rbind(AUC.TCR.PI.2023, AUC.CAR.2023)

AUC.fc <- data.frame(treatment =AUC.2023$treatment,
                             NFAT.AUC =AUC.2023$NFAT.AUC/AUC.2023$NFAT.AUC[5],
                             NFKB.AUC =AUC.2023$NFKB.AUC/AUC.2023$NFKB.AUC[5])

AUC.fc2 <- subset(AUC.fc, AUC.fc$treatment %in% c("CC","I0P0","I1.4P0_add" ,"I0P50_add","I1.4P50_add","CD19-1st","CD19-CD28" ,"CD19-WT-41BB"))

x <- reshape2::melt(AUC.fc2)

ggplot(data = AUC.fc2, aes(x = NFKB.AUC, y = NFAT.AUC, color = factor(treatment))) +
  themo_demo+
  geom_hline(aes(yintercept = 1.121), color="grey")+geom_vline(aes(xintercept = 1.231), color="grey")+
  geom_point(size=2, alpha=0.7)

# ggsave("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4B_dotplot_AUC.fc.pdf", width=4, height=2 )

AUC <- AUC.fc2

###@@@ Estimation

coTGs.DEG.TCR.PI$TCR.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[1]/AUC$NFAT.AUC[5])+ # NFATC
                                 2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[1]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$PIC.est<- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[3]/AUC$NFAT.AUC[5])+ # NFAT
                                2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[3]/AUC$NFKB.AUC[4])) # NFKB

dtx <- data.frame(symbol= coTGs.DEG.TCR.PI$Symbol, RefID = coTGs.DEG.TCR.PI$RefID,
                  PIC.effect.size = coTGs.DEG.TCR.PI$PIC.240_vs_NC.0-coTGs.DEG.TCR.PI$PIC.est,
                  TCR.effect.size = coTGs.DEG.TCR.PI$S.240_vs_NC.0-coTGs.DEG.TCR.PI$TCR.est)


###############################################################################################################################
#################                                           CAR                                               #################
###############################################################################################################################
load("/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/Fig4_RNA_DEseq2.by.group_CAR_2023.RData") # produced by Fig4_RNA_DESeq2.by.group.CAR
# load("//storage/labnas5/HuangWen5/2020_sequencing/CAR_J2nd26_seq_RNA/DESeq2.by.gropu.CAR.Jurkat.210623.RData")

colnames(DEG.all_actCAR_dataframe) <- c("symbol","RefID","zeta","WTBB","CD28","KRBB","TCR")

upregulated_CARS_2023 <- subset(DEG.all_actCAR_dataframe, DEG.all_actCAR_dataframe$zeta >1|DEG.all_actCAR_dataframe$WTBB>1 |DEG.all_actCAR_dataframe$CD28>1)
co_upregulated_CARS_2023 <- subset(upregulated_CARS_2023, upregulated_CARS_2023$symbol %in% coTGs$target)
length(unique(upregulated_CARS_2023$symbol))
length(unique(co_upregulated_CARS_2023$symbol))

DEG.all_actCAR_exp <- subset(DEG.all_actCAR_dataframe[,1:6], DEG.all_actCAR_dataframe$symbol %in% coTGs.DEG.TCR.PI$Symbol)
DEG.all_exp <- merge(coTGs.DEG.TCR.PI[,c(1:2,11,14)], DEG.all_actCAR_exp, by.x="RefID", by.y="RefID")

###@@@ Estimation 

coTGs.DEG.TCR.PI$zeta.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[6]/AUC$NFAT.AUC[5])+ # NFATC
                                   2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[6]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$CD28.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[7]/AUC$NFAT.AUC[5])+ # NFAT
                                  2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[7]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$WTBB.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[8]/AUC$NFAT.AUC[5])+ # NFAT
                                  2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[8]/AUC$NFKB.AUC[4])) # NFKB

DEG.all_for_effect.size <- merge(coTGs.DEG.TCR.PI[,c(1:2,15:19)], DEG.all_exp[,c(1:4,6:9)], by="RefID")


effect.size <- data.frame(symbol= DEG.all_for_effect.size$Symbol.x,RefID= DEG.all_for_effect.size$RefID,
                  PIC.effect.size = DEG.all_for_effect.size$PIC.240_vs_NC.0 - DEG.all_for_effect.size$PIC.est,
                  TCR.effect.size = DEG.all_for_effect.size$S.240_vs_NC.0 - DEG.all_for_effect.size$TCR.est,
                  
                  zeta.FC.exp.VS.est = DEG.all_for_effect.size$zeta - DEG.all_for_effect.size$zeta.est,
                  CD28.FC.exp.VS.est = DEG.all_for_effect.size$CD28 - DEG.all_for_effect.size$CD28.est,
                  WTBB.FC.exp.VS.est = DEG.all_for_effect.size$WTBB - DEG.all_for_effect.size$WTBB.est)

# write.csv(effect.size, "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/effect.size_2023.csv", row.names = F)
# write.csv(DEG.all_for_effect.size, "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/DEG.RNA.log2FC_2023.csv", row.names = F)
# 
# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4E_effect.size_noscle_2023.pdf", width=4, height=6 )

heatmap.effect.size <-
pheatmap(effect.size[,3:7], scale = "none", #cluster.all.annopeaks[,6:14], scale = "row", 
         treeheight_row = 40,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(200),
         cluster_cols = F, cluster_rows = T, clustering_method="ward.D2",
         cutree_rows = 2,
         show_rownames= F, show_colnames=T,
         cellwidth = 20 ,cellheight =2.5 ,main="effect size",
         fontsize=10)

effect.size$cluster <- cutree(heatmap.effect.size$tree_row, k = 2)
syb.1 <- subset(effect.size, cluster==1)
syb.2 <- subset(effect.size, cluster==2)
syb.2 <- subset(syb.2, !symbol %in% intersect(syb.1$symbol, syb.2$symbol)) 
length(unique(syb.1$symbol))
length(unique(syb.2$symbol))
# syb.1 <- syb.1[!duplicated(syb.1$Symbol),]
# syb.2 <- syb.2[!duplicated(syb.2$Symbol),]


# # "F5"      "IL2"     "PKM"     "NFKBID"  "GBP2"    "POU2AF1" "NR4A3"   "FXYD5"
# write.csv(syb.2, "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.effect.size.syb.2.csv", row.names = F)
# write.csv(syb.1, "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.effect.size.syb.1.csv", row.names = F)
# write.csv(unique(syb.2$symbol), "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.syb.2.csv", row.names = F)
# write.csv(unique(syb.1$symbol), "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.syb.1.csv", row.names = F)
# dev.off()

# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4E_RNA_effect.size.clustered_nonescale.pdf", width=4, height=6 )
# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4E_RNA_effect.size.clustered_rowscale.pdf", width=4, height=6 )

newOrder <- effect.size[heatmap.effect.size$tree_row$order,]
newOrder.DEG <- data.frame( matrix(NA, nrow(newOrder), ncol(DEG.all_for_effect.size)) )

for (gi in 1:nrow(newOrder)){
  newOrder.DEG[gi,] <- subset(DEG.all_for_effect.size, DEG.all_for_effect.size$Symbol.x == newOrder$symbol[gi])
}
colnames(newOrder.DEG) = colnames(DEG.all_for_effect.size)

pheatmap(newOrder.DEG[,c(9:14)], scale = "none", #cluster.all.annopeaks[,6:14], scale = "row", 
         treeheight_row = 40,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(200),
         cluster_cols = F, cluster_rows = F, #clustering_method="ward.D",
         cutree_rows = 2,
         gaps_row = nrow(syb.1),
         show_rownames= F, show_colnames=T,
         cellwidth = 20 ,cellheight =2.5 ,main="RNA (log2FC)",
         fontsize=10)
# dev.off()

library(corrplot)
library(RColorBrewer)

# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4C_effect.size_corrplot.pdf", width=7, height=7 )

M <- cor(effect.size[3:7])
p.mat <- cor.mtest(effect.size[3:7], conf.level = 0.95)

# display.brewer.all() 

## circle + colorful number
col3 = colorRampPalette(c('white', 'SteelBlue3', 'SteelBlue3','red'))

corrplot(M, order = "AOE", type = "upper", tl.pos = "d",p.mat=p.mat$p, insig = 'label_sig',
         is.corr = T, col.lim = c(-1,1), col = col3(100),
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2, pch.col = 'white')

corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",col = col3(100),
         diag = FALSE, tl.pos = "n", cl.pos = "n")

 # dev.off()
 
cor.TCR.PI <- cor.test(effect.size[,4], effect.size[,3]) #PI
cor.TCR.zeta <- cor.test(effect.size[,4], effect.size[,5]) #zeta
cor.TCR.CD28<- cor.test(effect.size[,4], effect.size[,6]) #CD28
cor.TCR.BB <- cor.test(effect.size[,4], effect.size[,7]) #BB
 
cor.list <- c("cor.TCR.PI", "cor.TCR.zeta", "cor.TCR.CD28","cor.TCR.BB")
p.list <- c("p.PI", "p.zeta", "p.CD28","p.WTBB" )
y.list <- c(3,5,6,7)

for (ti in 1:4){
  cor = get(cor.list[ti])
  plt <- data.frame(x=effect.size$TCR.effect.size, y= effect.size[,y.list[ti]])
  p <- ggplot(data = plt, aes(x = x , y =y )) + #PIC.effect.size)) + #
        geom_point(colour="navy", size=1, alpha=0.5)+
        themo_demo+ geom_abline(slope=1, intercept=0)+
        ylim(-6,6)+ xlim(-6,6)+
        labs(title= paste(p.list[ti],signif(cor$estimate,3), signif(cor$p.value,3), sep="_"), 
             x= "log2 FC TCR experiment / estimation",y= "log2 FC")+
         geom_hline(aes(yintercept = 0), color="grey")+ geom_vline(aes(xintercept = 0), color="grey")
  
  assign(paste("p",ti,sep=""), p)
}

###################
###  Figure 4D ####
###################   
# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4D_RNA_AUC_estimation.VS.exp_CARs.pdf", width =3,height =3)
ggpubr::ggarrange(p1,p3,p2,p4, nrow=2, ncol=2)
# dev.off()

###################
###  Figure 4F ####
###################
### /storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/GO/syb2_all.tw98r83_0
### /storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/GO/syb1_all.tzlt77le7

tmpdir <- "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/GO/"
GO.typ <- readxl::read_xlsx(paste(tmpdir,"syb2_all.tw98r83_0/metascape_result.xlsx",sep=""), sheet = "Enrichment") ##HiAUC_all.th2rpd5gk ##LOWAUC_all.t2q0hubfw
# GO.typ <- readxl::read_xlsx(paste(tmpdir,"syb1_all.tzlt77le7/metascape_result.xlsx",sep=""), sheet = "Enrichment") ##HiAUC_all.th2rpd5gk ##LOWAUC_all.t2q0hubfw

#
GO.typ <- GO.typ[grep("Summary",GO.typ$GroupID),]
GO.typ <- GO.typ[grep("GO Bio",GO.typ$Category),]

library(ggrepel)

for (nr in 1:nrow(GO.typ)){
  GO.typ$Genenumber[nr] = length(strsplit(GO.typ$Genes[nr], split=",")[[1]])
}

n=10
go_enrich_df<-data.frame(ID=GO.typ$GroupID[1:n],
                         Description=GO.typ$Description[1:n],
                         `-logP`= -(GO.typ$LogP)[1:n],
                         type=GO.typ$Category[1:n],
                         Genenumber = GO.typ$Genenumber[1:n])

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## add GO terms description names
labels <- go_enrich_df$Description
names(labels) = rev(1:nrow(go_enrich_df))


pa <-
# pb<-
  ggplot(data=go_enrich_df, aes(x=number, y=X.logP, fill=Genenumber)) +coord_flip() +
  geom_bar(stat="identity", width=0.5) +
  # scale_fill_manual(values = CPCOLS) +
  scale_fill_gradient(low="LightSkyBlue2",high="orangeRed")+
  theme_test() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  # labs(title = paste("top0.2_exp.VS.est", sep=""))
  labs(title = paste("eft.size.cluster", sep=""))

# pdf("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/Fig4F_eft.size.metascape_GO_exp.VS.est.pdf", width=20, height = 2)
ggarrange(pa,pb)
# dev.off()


####################### homer motif analysis ######################################################
# write.csv(unique(effect.size$symbol),"/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.eft.size.unique.symbols.csv", row.names=F)
eft.size.unique.symbols <- read.csv("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/2023.eft.size.unique.symbols.csv")

library(org.Hs.eg.db)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
columns <- c("tx_name", "gene_id")

PR <- getPromoters(txdb, upstream=5000, downstream=5000, by="transcript")
hg38.tss.5kb_bed <- data.frame(PR)

hg38.tss.5kb_bed <- subset(hg38.tss.5kb_bed, hg38.tss.5kb_bed$start > 0 & seqnames %in% c(paste("chr",1:22,sep=""), "chrX","chrY"))
hg38.tss.5kb_bed$peakID <- paste("peakID",1:nrow(hg38.tss.5kb_bed), sep="")
# # colnames(hg38.tss.5kb_bed)[1] <- "Chr"
pks <- GRanges(hg38.tss.5kb_bed)

peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                             tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
                             annoDb="org.Hs.eg.db")
pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)

syb58.effect.size.annotated <- subset(pks.annotated, pks.annotated$SYMBOL %in% eft.size.unique.symbols$x)
syb58.effect.size.annotated <- syb58.effect.size.annotated[!duplicated(syb58.effect.size.annotated$transcriptId),]
syb58.bed <- subset(hg38.tss.5kb_bed, hg38.tss.5kb_bed$peakID%in%  syb58.effect.size.annotated$peakID)
syb58.bed <- syb58.bed[order(syb58.bed$seqnames),]
syb58.bed <- data.frame(peakID = syb58.bed$peakID, Chr=syb58.bed$seqnames, Start=syb58.bed$start, End=syb58.bed$end, Strand=syb58.bed$strand)
# write.table(syb58.bed, "/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/homer/2023_syb58.bed.txt", col.names=F, row.names = F,quote=F,sep="\t")

## homer:annotatePeaks.pl 2023_syb58.bed.txt hg38 -m NFAT1_MA0152.1.motif -mask > motif_results_NFAT1.txt
## homer:annotatePeaks.pl 2023_syb58.bed.txt hg38 -m p65.motif -mask > motif_results_p65.txt
########################
##############
#### FigS8D ##
##############

# library(stringr)

filename <- dir("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/homer/Results/")
for (f in 1:length(filename)){
  ipt.homer <- read.csv(paste("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/homer/Results/", filename[f],sep=""), header = T)
  assign(strsplit(filename[f],split="[.]")[[1]][1], ipt.homer)
}


c1.p65 <- subset(motif_results_p65_2023, Gene.Name %in% syb.1$symbol)
c1.NFAT  <- subset(motif_results_NFAT1_2023, Gene.Name %in% syb.1$symbol)
c2.p65 <- subset(motif_results_p65_2023, Gene.Name %in% syb.2$symbol)
c2.NFAT  <- subset(motif_results_NFAT1_2023, Gene.Name %in% syb.2$symbol)

c1 <- merge(c1.p65, c1.NFAT, by="PeakID")
c2 <- merge(c2.p65, c2.NFAT, by="PeakID")


for (i in 1:2){
  ipt <- get(c("c1", "c2")[i])
  
  num.AB <- data.frame()
  mean.dt <- data.frame()
  for (j in 1:nrow(ipt)){
    bs.NFAT <- str_extract_all(ipt$NFAT.Homer.Distance.From.Peak[j], "\\d+")
    nka <- length(bs.NFAT[[1]])/3    
    bs.NFKB <- str_extract_all(ipt$p65.Homer.Distance.From.Peak[j], "\\d+")
    nkb <- length(bs.NFKB[[1]])/3
    num.AB[j,1:5] <- c(ipt$Gene.Name.x[j],ipt$PeakID[j], nka,nkb,c("c1", "c2")[i])
    
    bs.NFAT.new <- data.frame() 
    bs.NFKB.new <- data.frame() 
    for (k in 1:nka){ bs.NFAT.new[k,1] <- as.numeric(bs.NFAT[[1]][3*k-2]) }
    for (k in 1:nkb){ bs.NFKB.new[k,1] <- as.numeric(bs.NFKB[[1]][3*k-2]) }
    
    dt <- c()
    for (a in 1:nka){
      dti <- as.numeric(bs.NFAT.new$V1[a]-bs.NFKB.new$V1)
      dt <- c(dt, dti)}
    mean.dt[j,1:6] <- data.frame(ipt$Gene.Name.x[j], ipt$PeakID[j], mean(dt), mean(abs(dt)),min(abs(dt)), (c("c1", "c2")[i]))
  }
  assign(paste(c("c1", "c2")[i], ".mean.dt", sep=""),mean.dt)
  assign(paste("num.AB",c("c1", "c2")[i],sep=""), num.AB)
}

num.AB <- rbind(num.ABc1, num.ABc2)
num.AB$V3 <- as.numeric(num.AB$V3)
num.AB$V4 <- as.numeric(num.AB$V4)

colnames(c1.mean.dt)= colnames(c2.mean.dt)= c("Symbol","peakID","mean","mean.abs","min.abs","cluster")
cluster.mean.dt <- rbind(c1.mean.dt, c2.mean.dt)

# comparisons = compaired,
cluster.mean.dt$cluster <- factor(cluster.mean.dt$cluster , levels = c("c1","c2"))
min.abs.t <- t.test(c1.mean.dt$min.abs, c2.mean.dt$min.abs, alternative = "two.sided" , mu= 0, paired = FALSE, var.equal = FALSE,onf.level = 0.95)

ggplot(data = cluster.mean.dt) + geom_boxplot(aes(x = cluster, y = min.abs)) +themo_demo + ylim(-100,5000)+
  labs(title = paste(signif(min.abs.t$estimate[1],3), signif(min.abs.t$estimate[2],3),signif(min.abs.t$p.value, 3),sep="_"))

# ggsave("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/FigS8_Homer_AB_motif_min.distance.pdf")


colnames(num.AB) <- c("Symbol","PeakID","A","B","cluster")
# num.AB
plt.num.AB <- reshape2::melt(num.AB)

# compaired <- list(c("A","B"))
t.test(subset(plt.num.AB, cluster=="c1"& variable=="A")$value, subset(plt.num.AB, cluster=="c2"& variable=="A")$value)
Bnum.t <- t.test(subset(plt.num.AB, cluster=="c1"& variable=="B")$value, subset(plt.num.AB, cluster=="c2"& variable=="B")$value)

ggplot(data = plt.num.AB) + 
  geom_boxplot(aes(x = cluster, y = value, fill=variable)) +themo_demo +
  labs(title=paste(signif(Bnum.t$estimate[1],3),signif(Bnum.t$estimate[2],3), signif(Bnum.t$p.value,3), sep="_"),
       y="Motif count")

# ggsave("/storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/FigS8_AB_motif_num.pdf")

# ##########################                 ATAC-seq data              #################################################
# ##########################      combine heatmap in Fig5 with ATAC     #################################################
# ##########################            PTS bedtools multicov           #################################################

# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-30-1.last.bam rep2/renamed.bam/I-30-2.last.bam rep1/renamed.bam/P-30-1.last.bam rep2/renamed.bam/P-30-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-30-2.last.bam -bed merge_summits_peak_sep_rep_210428/bedops_merged_macs2.bed > ./merge_summits_peak_sep_rep_210428/pks30min_readcounts_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-120-1.last.bam rep2/renamed.bam/I-120-2.last.bam rep1/renamed.bam/P-120-1.last.bam rep2/renamed.bam/P-120-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-120-2.last.bam -bed merge_summits_peak_sep_rep_210428/bedops_merged_macs2.bed > ./merge_summits_peak_sep_rep_210428/pks120min_readcounts_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-240-1.last.bam rep2/renamed.bam/I-240-2.last.bam rep1/renamed.bam/P-240-1.last.bam rep2/renamed.bam/P-240-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-240-2.last.bam -bed merge_summits_peak_sep_rep_210428/bedops_merged_macs2.bed > ./merge_summits_peak_sep_rep_210428/pks240min_readcounts_bedtools_multicov.txt

X<- read.table("//storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/from last.bam/bedops_merged_macs2.bed")
colnames(X) <- c("Chr","Start","End","peakID","V")
pks <- GRanges(X[,1:4])
peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                             tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
                             annoDb="org.Hs.eg.db")
pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)

dir <- "//storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/from last.bam/"
tps <- c("pks30min_readcounts_bedtools_multicov.txt", "pks120min_readcounts_bedtools_multicov.txt", "pks240min_readcounts_bedtools_multicov.txt")

for (tpi in 1:3){
  pks_readcounts <- read.table(paste(dir, tps[tpi], sep=""), header = F)
  
  colnames(pks_readcounts) <- c("chr","start","end","peakID","V","NC.1","NC.2",
                                "I.1","I.2","P.1","P.2","PI.1","PI.2")
  
  Y <- data.frame(pks_readcounts[,c(1:4)],
                  NC= (pks_readcounts$NC.1+pks_readcounts$NC.2)/2,
                  I= (pks_readcounts$I.1+pks_readcounts$I.2)/2,
                  P= (pks_readcounts$P.1+pks_readcounts$P.2)/2,
                  PI= (pks_readcounts$PI.1+pks_readcounts$PI.2)/2)
  

  Y.annotated <- merge(Y, pks.annotated, by="peakID")
  
  Y.annotated.FC <- data.frame(Y.annotated[,1:4], I= Y.annotated$I/(Y.annotated$NC+0.1), 
                               P= Y.annotated$P/(Y.annotated$NC+0.1),  
                               PI= Y.annotated$PI/(Y.annotated$NC+0.1), Y.annotated[,c(9,17,19,20)] )
  
  assign(paste("FC.", tpi, sep=""), Y.annotated.FC)
}


FC <- data.frame( FC.1[,c(1:4,8:10)], FC.1[,5:7],  FC.2[,5:7],  FC.3[,5:7])
FC.coTGs <- subset(FC, FC$SYMBOL %in% coTGs$target)

AUC.I <- (1+ FC.coTGs$I)*30/2+( FC.coTGs$I+ FC.coTGs$I.1)*90/2+( FC.coTGs$I.1+ FC.coTGs$I.2)*120/2
AUC.P <- (1+ FC.coTGs$P)*30/2+( FC.coTGs$P+ FC.coTGs$P.1)*90/2+( FC.coTGs$P.1+ FC.coTGs$P.2)*120/2
AUC.PI <- (1+ FC.coTGs$PI)*30/2+( FC.coTGs$PI+ FC.coTGs$PI.1)*90/2+( FC.coTGs$PI.1+ FC.coTGs$PI.2)*120/2

AUC.FC.coTGs <- data.frame( FC.coTGs[,1:7], AUC.I,  AUC.P,  AUC.PI)
symbol.for.AUC.FC.coTGs <- data.frame(symbol= unique(AUC.FC.coTGs$SYMBOL))

maxpk.AUC.FC.coTGs<- data.frame()
for (syb in 1:nrow(symbol.for.AUC.FC.coTGs)){
  x <- subset(AUC.FC.coTGs, SYMBOL == symbol.for.AUC.FC.coTGs[syb,1])
  maxpk.AUC.FC.coTGs <- rbind(maxpk.AUC.FC.coTGs, x[which.max(x$AUC.PI),])
}


COMP <- list(c(" AUC.PI"," AUC.I"),c(" AUC.PI"," AUC.P"))
### all coTGs
ptdata <- reshape2::melt( maxpk.AUC.FC.coTGs[,c(8:10)],)
b1 <-
  ggplot(ptdata, aes(y=as.numeric(value), x=variable, fill=variable))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  geom_signif(comparisons =COMP ,step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo+labs(title =paste("AUC of ",nrow(maxpk.AUC.FC.coTGs), " coTGs", sep=""))+
  geom_hline(aes(yintercept = 240), color="grey")
t.test(maxpk.AUC.FC.coTGs[,8], maxpk.AUC.FC.coTGs[,10])
t.test(maxpk.AUC.FC.coTGs[,9], maxpk.AUC.FC.coTGs[,10])

### PI up-regulated coTGs
pt <- subset( maxpk.AUC.FC.coTGs, SYMBOL%in% coDEG.PI$Symbol)
pt2 <- reshape2::melt(pt[,c(8:10)],)
b2 <-
  ggplot(pt2, aes(y=as.numeric(value), x=variable, fill=variable))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  geom_signif(comparisons =COMP ,step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo+labs(title =paste("AUC of ",nrow(pt)," coTGs up-regulated under PMA&Iono",sep=""))+
  geom_hline(aes(yintercept = 240), color="grey")
t.test(pt$AUC.PI, pt$AUC.P)
t.test(pt$AUC.PI, pt$AUC.I)
t.test(pt$AUC.PI, maxpk.AUC.FC.coTGs$AUC.PI)


## 2 group coTGs
maxpk.AUC.FC.syb1 <- subset( maxpk.AUC.FC.coTGs, SYMBOL %in% syb.1$symbol)
maxpk.AUC.FC.syb1$group= "group1"
maxpk.AUC.FC.syb2 <- subset( maxpk.AUC.FC.coTGs, SYMBOL %in% syb.2$symbol)
maxpk.AUC.FC.syb2$group= "group2"

maxpk.AUC.FC.syball <- rbind( maxpk.AUC.FC.syb2,  maxpk.AUC.FC.syb1)

pt2<- reshape2::melt( maxpk.AUC.FC.syball[,c(8:11)],)
b3 <-
  ggplot(pt2, aes(y=as.numeric(value), x=variable, fill=group))+
  geom_boxplot(notch = F, notchwidth = 0.5)+  
  geom_signif(comparisons =COMP ,step_increase = 0.1, map_signif_level = F,test = t.test)+
  # geom_signif(comparisons =list(c("group1","group2")), step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo+labs(title =paste("AUC of ",nrow(maxpk.AUC.FC.syball)," coTGs in 2 groups", sep=""))+
  geom_hline(aes(yintercept = 240), color="grey")


t.test(subset( maxpk.AUC.FC.syball, group=="group1")[,10], subset( maxpk.AUC.FC.syball, group=="group2")[,10])
t.test(subset( maxpk.AUC.FC.syball, group=="group1")[,9], subset( maxpk.AUC.FC.syball, group=="group2")[,9])
t.test(subset( maxpk.AUC.FC.syball, group=="group1")[,8], subset( maxpk.AUC.FC.syball, group=="group2")[,8])

# pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/FigS8_HIJ_maxpk.AUC.FC.pdf",width=10,height=4)
ggpubr::ggarrange(b1,b2,b3, nrow=1, ncol=3)
# 
# dev.off()

gs <- c("SH2D2A","CD69","CD27")
for (i in 1:3){
  gi<-gs[i]
  
  
  sp <- subset(FC.coTGs, SYMBOL== gi)[,c(1,8:16)]
  sp <- subset(sp, sp$peakID %in% maxpk.AUC.FC.coTGs$peakID)
  for(pkid in 1:nrow(sp)){
    sp2 <- subset(sp, peakID == sp$peakID[pkid])
    label = data.frame(time=as.numeric(rep(c(0,30,120,240),3)), group=c(rep("Iono",4),rep("PMA",4),rep("PMA & Iono",4)))
    label$value <- as.numeric(c(1,sp2$I,sp2$I.1,sp2$I.1, 
                                1,sp2$P,sp2$P.1,sp2$P.2,
                                1,sp2$PI,sp2$PI.1,sp2$PI.2))
    # dt.line = reshape2::melt(label, id=c("time","group"),  measure.vars=c("value"))
    
    p<-
      ggplot(label, aes(x= time, y = value, color=group)) +
      geom_line(size = 1) +geom_point(size=3)+ themo_demo +
      labs(title=paste(gi,sep=""),y="Fold change",x="Time(min)")+ 
      scale_color_manual(values = c("SteelBlue3", "#FF8247","#7A378B"))+
      scale_x_continuous(breaks = c(0, 30, 120, 240))#+ylim(-5,5)
    
    assign(paste("p",i,pkid,sep="") ,p)
    
  }
  
}
# 
pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig4/FigS8_G_maxpk.examples.FC.pdf",width=10,height=4)
ggpubr::ggarrange(p11,p21,p31, nrow=1, ncol=3)
dev.off()