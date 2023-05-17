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

###### load all required data before RNA-seq data analysis #### 20210923-report ####################################################
load("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/CD4_TF_TG_analysis_Step1_2023.RData")

############################################ RNA-seq data (DEG) #####################################################################
load("/storage/disk2/HuangW/htseq_count_results/DEseq2.by.group.210628_Ab_PI.RData") ## produced by "U:/HuangW/htseq_count_results/DEseq2.by.group.210628_Ab_PI.R"

trd <- log2(2)
UP_DEG.TCR <- subset(DEG.all_act.up_dataframe, S.240_vs_NC.0 > trd) ##(DEG.all.up_dataframe_exon, log2FC_S.240 > trd)
UP_DEG.PI <- subset(DEG.all_act.up_dataframe, PIC.240_vs_NC.0> trd) ##(DEG.all.up_dataframe_exon, log2FC_PIC.240> trd)
# UP_DEG.PI_P_I <- subset(DEG.all_act.up_dataframe, PIC.240_vs_NC.0> trd & I.240_vs_NC.0> 0 & P.240_vs_NC.0> 0) ##DEG.all.up_dataframe_exon
UP_DEGs <- unique(rbind(UP_DEG.TCR, UP_DEG.PI))

y1 <- calculate.overlap(x= list(TCR = unique(UP_DEG.TCR$Symbol),
                                PI = unique(UP_DEG.PI$Symbol)))

coDEG.TCR <- subset(UP_DEG.TCR, UP_DEG.TCR$Symbol %in% coTGs$target)
coDEG.PI <- subset(UP_DEG.PI, UP_DEG.PI$Symbol %in% coTGs$target)

y2 <- calculate.overlap(x= list(TCR = unique(coDEG.TCR$Symbol),
                                PI = unique(coDEG.PI$Symbol)))


UP_DEG.TCR.PI <- subset(UP_DEGs, UP_DEGs$Symbol %in% y1[["a3"]]) # 705 row, 669 gene
# coTGs.DEG.TCR.PI1 <- subset(UP_DEG.TCR.PI, UP_DEG.TCR.PI$Symbol %in% coTGs$target) #232 row, 218 gene
# coTGs.DEG.TCR.PI2 <- subset(UP_DEG.TCR.PI, UP_DEG.TCR.PI$Symbol %in% coTG.dorothea$target) #232row,218gene
# coTGs.DEG.TCR.PI <- unique(rbind(coTGs.DEG.TCR.PI1, coTGs.DEG.TCR.PI2))
coTGs.DEG.TCR.PI <- subset(UP_DEG.TCR.PI, UP_DEG.TCR.PI$Symbol %in% coTGs$target) #104 row, 96 gene
coTGs.DEG.TCR.PI <- coTGs.DEG.TCR.PI[!duplicated(coTGs.DEG.TCR.PI$Symbol),]


# write.csv(coTGs.DEG.TCR.PI, paste(tmpdir,"/combinatory_effect/coTGs.DEG.TCR.PI.csv",sep=""))



######## from D:/LabNAS_Backup/HuangWen/Movie/Dragonfly/2020/PI-matrix-scdata-from-20200813-0814/11report_TCR_PI_T4h_analysis_HW.R
# AUC.nobasal <- read.csv("/storage/disk2/HuangW/HW_PI_timepoint_ATAC/chromVAR_code/11report_AUC.mean.unit.nobasal.csv")
# AUC.fc <- read.csv("/storage/disk2/HuangW/HW_PI_timepoint_ATAC/chromVAR_code/11report_AUC.mean.unit.fc.csv")
# AUC.fc <- read.csv("/storage/disk2/HuangW/HW_PI_timepoint_ATAC/chromVAR_code/11report_new_AUC.mean.unit.fc.csv")

AUC.TCR.PI.2023 <- read.csv("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/Fig4.AUC_mean_TCR_PI.csv")
AUC.CAR.2023 <- read.csv("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/Fig4.AUC_mean_CAR.csv")
AUC.2023 <- rbind(AUC.TCR.PI.2023, AUC.CAR.2023)

AUC.fc <- data.frame(treatment =AUC.2023$treatment,
                             NFAT.AUC =AUC.2023$NFAT.AUC/AUC.2023$NFAT.AUC[5],
                             NFKB.AUC =AUC.2023$NFKB.AUC/AUC.2023$NFKB.AUC[5])

AUC.fc2 <- subset(AUC.fc, AUC.fc$treatment %in% c("CC","I0P0","I1.4P0_add" ,"I0P50_add","I1.4P50_add","CD19-1st","CD19-CD28" ,"CD19-WT-41BB"))
# AUC.fc2 <- subset(AUC.fc, AUC.fc$treatment %in% c("CC","I0P0","I1P0" ,"I0P50","I1P50","CD19-1st","CD19-CD28" ,"CD19-WT-41BB"))

x <- reshape2::melt(AUC.fc)

ggplot(data = AUC.fc2, aes(x = NFKB.AUC, y = NFAT.AUC, color = factor(treatment))) +
  themo_demo+
  geom_hline(aes(yintercept = 1.121), color="grey")+geom_vline(aes(xintercept = 1.231), color="grey")+
  geom_point(size=2, alpha=0.7)

# ggsave("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/11reprot_dotplot_AUC.fc.pdf", width=4, height=2 )

AUC <- AUC.fc2
######## from D:/LabNAS_Backup/HuangWen/Movie/Dragonfly/2020/PI-matrix-scdata-from-20200813-0814/11report_TCR_PI_T4h_analysis_HW.R
###@@@ Estimation

coTGs.DEG.TCR.PI$TCR.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[1]/AUC$NFAT.AUC[3])+ # NFATC
                                 2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[1]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$PIC.est<- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[5]/AUC$NFAT.AUC[3])+ # NFAT
                                2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[5]/AUC$NFKB.AUC[4])) # NFKB

dtx <- data.frame(symbol= coTGs.DEG.TCR.PI$Symbol,
                  PIC.effect.size = coTGs.DEG.TCR.PI$PIC.240_vs_NC.0-coTGs.DEG.TCR.PI$PIC.est,
                  TCR.effect.size = coTGs.DEG.TCR.PI$S.240_vs_NC.0-coTGs.DEG.TCR.PI$TCR.est)


# # ################################################# CAR ##############################################################
# # ### run "W:\HuangWen5\2020_sequencing\CAR_J2nd26_seq_RNA\DESeq2.by.gropu.CAR.Jurkat.210623.R"
load("//storage/labnas5/HuangWen5/2020_sequencing/CAR_J2nd26_seq_RNA/DESeq2.by.gropu.CAR.Jurkat.210623.RData")
colnames(DEG.all_actCAR_dataframe) <- c("symbol","zeta","WTBB","CD28","KRBB","TCR")

upregulated_CARS_2023 <- subset(DEG.all_actCAR_dataframe[,1:4], DEG.all_actCAR_dataframe$zeta >1|DEG.all_actCAR_dataframe$WTBB>1 |DEG.all_actCAR_dataframe$CD28>1)
co_upregulated_CARS_2023 <- subset(upregulated_CARS_2023, upregulated_CARS_2023$symbol %in% coTGs$target)

DEG.all_actCAR_exp <- subset(DEG.all_actCAR_dataframe[,1:4], DEG.all_actCAR_dataframe$symbol %in% coTGs.DEG.TCR.PI$Symbol)
DEG.all_exp <- merge(coTGs.DEG.TCR.PI[,c(2,11,14)], DEG.all_actCAR_exp, by.x="Symbol", by.y="symbol")
###@@@ Estimation 

coTGs.DEG.TCR.PI$zeta.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[6]/AUC$NFAT.AUC[3])+ # NFATC
                                   2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[6]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$CD28.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[7]/AUC$NFAT.AUC[3])+ # NFAT
                                  2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[7]/AUC$NFKB.AUC[4])) # NFKB

coTGs.DEG.TCR.PI$WTBB.est <- log2(2^(coTGs.DEG.TCR.PI$I.240_vs_NC.0)*(AUC$NFAT.AUC[8]/AUC$NFAT.AUC[3])+ # NFAT
                                  2^(coTGs.DEG.TCR.PI$P.240_vs_NC.0)*(AUC$NFKB.AUC[8]/AUC$NFKB.AUC[4])) # NFKB

DEG.all_for_effect.size <- merge(coTGs.DEG.TCR.PI[,c(2,15:19)], DEG.all_exp, by="Symbol")


effect.size <- data.frame(symbol= DEG.all_for_effect.size$Symbol,
                  PIC.effect.size = DEG.all_for_effect.size$PIC.240_vs_NC.0 - DEG.all_for_effect.size$PIC.est,
                  TCR.effect.size = DEG.all_for_effect.size$S.240_vs_NC.0 - DEG.all_for_effect.size$TCR.est,
                  
                  zeta.FC.exp.VS.est = DEG.all_for_effect.size$zeta - DEG.all_for_effect.size$zeta.est,
                  CD28.FC.exp.VS.est = DEG.all_for_effect.size$CD28 - DEG.all_for_effect.size$CD28.est,
                  WTBB.FC.exp.VS.est = DEG.all_for_effect.size$WTBB - DEG.all_for_effect.size$WTBB.est)

# write.csv(effect.size, "/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.effect.size.csv", row.names = F)
# write.csv(DEG.all_for_effect.size, "/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.DEG.RNA.log2FC.csv", row.names = F)
# pdf("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.effect.size_noscle.pdf", width=4, height=6 )

heatmap.effect.size <-
pheatmap(effect.size[,2:6], scale = "none", #cluster.all.annopeaks[,6:14], scale = "row", 
         treeheight_row = 40,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(200),
         cluster_cols = F, cluster_rows = T, clustering_method="ward.D",
         cutree_rows = 2,
         show_rownames= F, show_colnames=T,
         cellwidth = 20 ,cellheight =2.5 ,main="effect size",
         fontsize=10)

effect.size$cluster <- cutree(heatmap.effect.size$tree_row, k = 2)
syb.1 <- subset(effect.size, cluster==1)
syb.2 <- subset(effect.size, cluster==2)
syb.2 <- subset(syb.2, !symbol %in% intersect(syb.1$symbol, syb.2$symbol))
# write.csv(syb.2, "/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.effect.size.syb.2.csv", row.names = F)
# write.csv(syb.1, "/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.effect.size.syb.1.csv", row.names = F)

# pdf("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/2023.DEG.all_for_effect.size.clustered_rowscale.pdf", width=4, height=6 )

newOrder <- effect.size[heatmap.effect.size$tree_row$order,]
newOrder.DEG <- data.frame( matrix(NA, nrow(newOrder), ncol(DEG.all_for_effect.size)) )

for (gi in 1:nrow(newOrder)){
  newOrder.DEG[gi,] <- subset(DEG.all_for_effect.size, DEG.all_for_effect.size$Symbol == newOrder$symbol[gi])
}
colnames(newOrder.DEG) = colnames(DEG.all_for_effect.size)

pheatmap(newOrder.DEG[,c(7,8,9,11,10)], scale = "none", #cluster.all.annopeaks[,6:14], scale = "row", 
         treeheight_row = 40,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(200),
         cluster_cols = F, cluster_rows = F, #clustering_method="ward.D",
         cutree_rows = 2,
         gaps_row = nrow(syb.2),
         show_rownames= F, show_colnames=T,
         cellwidth = 20 ,cellheight =2.5 ,main="RNA (log2FC)",
         fontsize=10)
dev.off()

library(corrplot)
library(RColorBrewer)

# pdf("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/effect.size_corrplot.pdf", width=4, height=6 )

 M <- cor(effect.size[2:6])
 p.mat <- cor.mtest(effect.size[2:6], conf.level = 0.95)
 
 
 # display.brewer.all() 
 
 ## circle + colorful number
 col3 = colorRampPalette(c('white', 'SteelBlue3', 'SteelBlue3','red'))
 
 corrplot(M, order = "AOE", type = "upper", tl.pos = "d",p.mat=p.mat$p, insig = 'label_sig',
          is.corr = T, col.lim = c(-1,1), col = col3(100),
          sig.level = c(0.001, 0.01, 0.05), pch.cex = 2, pch.col = 'white')
 
 corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",col = col3(100),
          diag = FALSE, tl.pos = "n", cl.pos = "n")
 
 # dev.off()
 
 cor.test(effect.size[,3], effect.size[,2])
 cor.test(effect.size[,3], effect.size[,4])
 cor.test(effect.size[,3], effect.size[,5])
 cor.test(effect.size[,3], effect.size[,6])
 
 p.zeta <-
   # p.CD28 <-
   # p.WTBB <-
   # p.PI <-
   ggplot(data = effect.size, aes(x = TCR.effect.size , y = zeta.FC.exp.VS.est)) +
   themo_demo+ geom_abline(slope=1, intercept=0)+
   ylim(-6,6)+ xlim(-6,6)+
   labs(title= paste(x$estimate, x$p.value, sep=" "), x= "log2 FC TCR experiment / estimation",y= "log2 FC")+
   geom_hline(aes(yintercept = 0), color="grey")+ geom_vline(aes(xintercept = 0), color="grey")+
   geom_point(colour="navy", size=1, alpha=0.5)
 
# pdf("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/RNA_AUC_estimation.VS.exp_CARs.pdf", width =8,height = 2)
  ggpubr::ggarrange(p.PI,p.zeta,p.CD28,p.WTBB, nrow=1, ncol=4)
# dev.off()

  

#  U:\HuangW\T_Cell_Dynamic_paper\metadata\eft.size.cluster2_all.tab84u8_b
#  U:\HuangW\T_Cell_Dynamic_paper\metadata\eft.size.cluster1_all.tvkhp8rzq
  
tmpdir <- "//storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/"
GO.typ <- readxl::read_xlsx(paste(tmpdir,"eft.size.cluster1_all.tvkhp8rzq/metascape_result.xlsx",sep=""), sheet = "Enrichment") ##HiAUC_all.th2rpd5gk ##LOWAUC_all.t2q0hubfw
# GO.typ <- readxl::read_xlsx(paste(tmpdir,"Low_AUC_all.tw2ku3nbn/metascape_result.xlsx",sep=""), sheet = "Enrichment") ##HiAUC_all.th2rpd5gk ##LOWAUC_all.t2q0hubfw

#
GO.typ <- GO.typ[grep("Summary",GO.typ$GroupID),]
GO.typ <- GO.typ[grep("GO Bio",GO.typ$Category),]

# OP_gene <- data.frame(strsplit(GO.typ$Symbols[1],split="[,]"))
# OP_gene_yyyy <- subset(yyyy, yyyy$symbol %in% OP_gene[,1])

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


pb <- ggplot(data=go_enrich_df, aes(x=number, y=X.logP, fill=Genenumber)) +coord_flip() +
  geom_bar(stat="identity", width=0.5) +
  # scale_fill_manual(values = CPCOLS) +
  scale_fill_gradient(low="LightSkyBlue2",high="orangeRed")+
  theme_test() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  # labs(title = paste("top0.2_exp.VS.est", sep=""))
  labs(title = paste("eft.size.cluster1", sep=""))

# pdf(paste(tmpdir,"eft.size.cluster1_metascape_GO_exp.VS.est.pdf"), width=8, height = 2)
# ggarrange(pa,pb)
pb
# dev.off()

# pdf(paste(tmpdir,"RNA_Low_AUC_all.tw2ku3nbn_metascape_GO_exp.VS.est.pdf"), width=8, height =8)
# # ggarrange(pa,pb)
# pb
# dev.off()

####################### homer motif analysis ######################################################
library(stringr)

filename <- dir("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/homer_cluster12/2023")
for (f in 1:length(filename)){
  ipt.homer <- read.csv(paste("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/homer_cluster12/2023/", filename[f],sep=""), header = T)
  assign(strsplit(filename[f],split="[.]")[[1]][1], ipt.homer)
}

# c1 <- merge(cluster1_annoNFATC1, cluster1_annop65, by="PeakID")
# c2 <- merge(cluster2_annoNFATC1, cluster2_annop65, by="PeakID")
c1 <- merge(`2023cluster1_annoNFAT1`, cluster1_annop65, by="PeakID")
c2 <- merge(`2023cluster2_annoNFAT1`, cluster2_annop65, by="PeakID")

num.AB <- data.frame()
for (i in 1:2){
  ipt <- get(c("c1", "c2")[i])
  
  mean.dt <- data.frame()
  for (j in 1:nrow(ipt)){
    bs.NFAT <- str_extract_all(ipt$NFAT.Homer.Distance.From.Peak[j], "\\d+")
    nka <- length(bs.NFAT[[1]])/3    
    bs.NFKB <- str_extract_all(ipt$p65.Homer.Distance.From.Peak[j], "\\d+")
    nkb <- length(bs.NFKB[[1]])/3
    num.AB[j,1:4] <- data.frame(ipt$Gene.Name.x[j],nka,nkb,c("c1", "c2")[i])
    
    bs.NFAT.new <- data.frame() 
    bs.NFKB.new <- data.frame() 
    for (k in 1:nka){ bs.NFAT.new[k,1] <- as.numeric(bs.NFAT[[1]][3*k-2]) }
    for (k in 1:nkb){ bs.NFKB.new[k,1] <- as.numeric(bs.NFKB[[1]][3*k-2]) }
    
    dt <- c()
    for (a in 1:nka){
      dti <- as.numeric(bs.NFAT.new$V1[a]-bs.NFKB.new$V1)
      dt <- c(dt, dti)}
    mean.dt[j,1:5] <- data.frame(ipt$Gene.Name.x[j], mean(dt), mean(abs(dt)),min(abs(dt)), (c("c1", "c2")[i]))
  }
  assign(paste(c("c1", "c2")[i], ".mean.dt", sep=""),mean.dt)
}

colnames(c1.mean.dt)= colnames(c2.mean.dt)= c("Symbol","mean","mean.abs","min.abs","cluster")
cluster.mean.dt <- rbind(c1.mean.dt, c2.mean.dt)

# comparisons = compaired,
cluster.mean.dt$cluster <- factor(cluster.mean.dt$cluster , levels = c("c2","c1"))
ggplot(data = cluster.mean.dt) + geom_boxplot(aes(x = cluster, y = min.abs)) +themo_demo + ylim(-100,5000)
t.test(c1.mean.dt$min.abs, c2.mean.dt$min.abs, alternative = "two.sided" , mu= 0, paired = FALSE, var.equal = FALSE,onf.level = 0.95)

ggplot(cluster.mean.dt, aes(x=mean.abs, fill=cluster)) +themo_demo+
  geom_density(alpha=.25)

# ggsave("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/homer_cluster12/2023_Home_AB_motif_min.distance.pdf")

colnames(num.AB) <- c("Symbol","A","B","cluster")
plt.num.AB <- reshape2::melt(num.AB)

# compaired <- list(c("A","B"))
ggplot(data = plt.num.AB) + 
  geom_boxplot(aes(x = cluster, y = value, fill=variable)) +themo_demo 
# geom_signif(comparisons = compaired,step_increase = 0.1, map_signif_level = F,test = t.test)
t.test(subset(plt.num.AB, cluster=="c1"& variable=="A")$value, subset(plt.num.AB, cluster=="c2"& variable=="A")$value)
t.test(subset(plt.num.AB, cluster=="c1"& variable=="B")$value, subset(plt.num.AB, cluster=="c2"& variable=="B")$value)

# ggsave("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/homer_cluster12/2022_Home_AB_motif_num.pdf")


# ############################################ ATAC-seq data ################################################################
# ##########################        combine heatmap in Fig5 with ATAC
# ###########################################################################################################################

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# tn5.hit.counts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/pks_readcounts_bedtools_multicov.txt", header = T)
# #pks_tn5 <- write.table(tn5.hit.counts[,c(2:4)], "/storage/labnas5/HuangWen5/linwei/linux/pks_in_counts_hw.bed", sep = "\t",
#                        row.names = F, col.names = F, quote=F)
# pks <- GRanges(tn5.hit.counts[,c(2:5,1)])
# peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
#                              tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
#                              annoDb="org.Hs.eg.db")
# pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)
# tn5.hit.syb2 <- subset(pks.annotated, SYMBOL %in% syb.2$symbol)
# tn5.hit.syb2.counts <- subset(tn5.hit.counts, Geneid %in% tn5.hit.syb2$Geneid)
# tn5.hit.syb2.counts$group= "group1"
# tn5.hit.syb1 <- subset(pks.annotated, SYMBOL %in% syb.1$symbol)
# tn5.hit.syb1.counts <- subset(tn5.hit.counts, Geneid %in% tn5.hit.syb1$Geneid)
# tn5.hit.syb1.counts$group= "group2"
#
# tn5.hit.syball.counts <- rbind(tn5.hit.syb2.counts, tn5.hit.syb1.counts)
#
# tn5.hit.syball.meancounts <- data.frame(tn5.hit.syball.counts[,c(1:6,15)],
#                                       NC=(tn5.hit.syball.counts$NC.0.1.last.bam+tn5.hit.syball.counts$NC.0.2.last.bam)/2,
#                                       I=(tn5.hit.syball.counts$I.240.1.last.bam+tn5.hit.syball.counts$I.240.2.last.bam)/2,
#                                       P=(tn5.hit.syball.counts$P.240.1.last.bam+tn5.hit.syball.counts$P.240.2.last.bam)/2,
#                                       PI=(tn5.hit.syball.counts$PI.240.1.last.bam+tn5.hit.syball.counts$PI.240.2.last.bam)/2)
#
# tn5.hit.syball.FC <- data.frame(tn5.hit.syball.meancounts[,1:7],
#                                       I= tn5.hit.syball.meancounts$I/(tn5.hit.syball.meancounts$NC+1),
#                                       P= tn5.hit.syball.meancounts$P/(tn5.hit.syball.meancounts$NC+1),
#                                       PI= tn5.hit.syball.meancounts$PI/(tn5.hit.syball.meancounts$NC+1))
#
# nu <- nrow(tn5.hit.syb2.counts)
# pheatmap(tn5.hit.syball.FC[,8:10], scale = "row", #cluster.all.annopeaks[,6:14], scale = "row",
#          # treeheight_row = 40,
#          color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(200),
#          cluster_cols = F, cluster_rows = F, #clustering_method="ward.D",
#          gaps_row=nu,
#          show_rownames= F, show_colnames=T,
#          cellwidth = 20 ,cellheight =0.5 ,main="row scaled FC of peak counts",
#          fontsize=10)
#
# ptdata <- reshape2::melt(tn5.hit.syball.FC[,7:10])
# ggplot(ptdata, aes(y=value, x=variable, fill=group))+
#   geom_boxplot(notch = T, notchwidth = 0.5)+
#   themo_demo+labs(title ="peak-count-FC of group1&2 coTGs (Tn5 points)")+
#   geom_hline(aes(yintercept = 1), color="grey")+
#   ylim(0,5)
#
# t.test(subset(tn5.hit.syball.FC, group=="group1")[,10], subset(tn5.hit.syball.FC, group=="group2")[,10])
# t.test(subset(tn5.hit.syball.FC, group=="group1")[,9], subset(tn5.hit.syball.FC, group=="group2")[,9])
# t.test(subset(tn5.hit.syball.FC, group=="group1")[,8], subset(tn5.hit.syball.FC, group=="group2")[,8])
#
# ########## tn5 counts and coDEGs
#
# tn5.coTGs.anno <- subset(pks.annotated, pks.annotated$SYMBOL %in% coTGs$target)
# tn5.coTGs.counts <- subset(tn5.hit.counts, tn5.hit.counts$Geneid %in% tn5.coTGs.anno$Geneid)
#
# tn5.coTGs.meancounts <- data.frame(tn5.coTGs.counts[,c(1:6)],
#                                         NC=(tn5.coTGs.counts$NC.0.1.last.bam+tn5.coTGs.counts$NC.0.2.last.bam)/2,
#                                         I=(tn5.coTGs.counts$I.240.1.last.bam+tn5.coTGs.counts$I.240.2.last.bam)/2,
#                                         P=(tn5.coTGs.counts$P.240.1.last.bam+tn5.coTGs.counts$P.240.2.last.bam)/2,
#                                         PI=(tn5.coTGs.counts$PI.240.1.last.bam+tn5.coTGs.counts$PI.240.2.last.bam)/2)
#
# tn5.coTGs.FC <- data.frame(tn5.coTGs.meancounts[,1:6],
#                                 I= tn5.coTGs.meancounts$I/(tn5.coTGs.meancounts$NC+0.001),
#                                 P= tn5.coTGs.meancounts$P/(tn5.coTGs.meancounts$NC+0.001),
#                                 PI= tn5.coTGs.meancounts$PI/(tn5.coTGs.meancounts$NC+0.001))
#
#
# ptdata <- reshape2::melt(tn5.coTGs.FC[,7:9])
# ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
#   geom_boxplot(notch = T, notchwidth = 0.5)+
#   themo_demo+labs(title ="Fold change of peak counts")+
#   geom_hline(aes(yintercept = 1), color="grey")+
#   ylim(0,5)
#
# ggplot(ptdata, aes(x=value, color=variable)) +
#   geom_density()+
#   xlim(0,4)+ labs(title=paste("Peak-count-FC distribution of all",
#                               "NFAT1&RELA co-targtes under 3 stimulations (Tn5 points)",sep="\n"))+
#   themo_demo
#
# ks.test(pks.coTGs.FC[,5], pks.coTGs.FC[,6])
# ks.test(pks.coTGs.FC[,5], pks.coTGs.FC[,7])
# ks.test(pks.coTGs.FC[,6], pks.coTGs.FC[,7])
# # ---------------------------------------------------------------------------
# # EdgeR
# # differential expressed gene analysis by edgeR
# # ---------------------------------------------------------------------------
# 
# library(ggplot2)
# library(reshape2)
# library(ggbeeswarm)
# library("edgeR")
# library(dplyr)
# library(grid)
# library(biomaRt)
# library(org.Hs.eg.db)
# 
# dir_save = "/storage/labnas5/HuangWen5/linwei/"
# setwd(dir_save)
# mycounts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/counts.txt", header = T)
# 
# rownames(mycounts)<-mycounts[,1]
# matrix <- mycounts[,7:14]
# colnames()
# #分组
# group_list =c(rep("I",2),rep("NC",2),rep("P",2),rep("PI",2))
# 
# counts
# #差异分析
# y <- DGEList(counts=matrix, group = group_list)
# y
# y <- calcNormFactors(y)
# y <- estimateCommonDisp(y)
# y <- estimateTagwiseDisp(y)
# 
# 
# ########################寻找差异基因##############################
# 
# #设计样本matrix
# design <- model.matrix(~0+group, data=y$samples)
# colnames(design) <- levels(y$samples$group)
# 
# fit <- glmQLFit(y, design)
# 
# #一对一比较
# qlf.I_NC  <- glmQLFTest(fit, contrast= c(1,-1,0,0))
# qlf.P_NC <- glmQLFTest(fit, contrast= c(0,-1,1,0))
# qlf.PI_NC <- glmQLFTest(fit, contrast=c(0,-1,0,1))
# # qlf.PI_I  <- glmQLFTest(fit, contrast= c(-1,0,0,1))
# # qlf.PI_P  <- glmQLFTest(fit, contrast= c(0,0,-1,1))
# # qlf.P_I  <- glmQLFTest(fit, contrast=c(-1,0,1,0))
# 
# #同时比较
# my.contrasts <- makeContrasts(P.NC =P-NC, I.NC =I-NC, PI.NC =PI-NC,
#                               # PI.I =PI-I, PI.P =PI-P, P.I =P-I,
#                               levels=design)
# 
# qlf<- glmQLFTest(fit, contrast=my.contrasts)
# #Example: qlf.S_10minvsNC <- glmQLFTest(fit, contrast=my.contrasts[,"S_10minvsNC"])
# 
# 
# #不筛选差异基�?
# tTag <- topTags(qlf, n=nrow(y))
# DEG <- tTag$table
# DEG.pks <- subset(mycounts, mycounts$Geneid%in% rownames(DEG))[,1:5]
# pks <- GRanges(DEG.pks[,c(2:5,1)])
# peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
#                              tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
#                              annoDb="org.Hs.eg.db")
# pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)
# 
# 
# tn5.hit.syb2 <- subset(pks.annotated, SYMBOL %in% syb.2$symbol)
# tn5.hit.syb2.DEG <- subset(DEG, rownames(DEG) %in% tn5.hit.syb2$Geneid)
# tn5.hit.syb2.DEG$group <- "group1"
# 
# tn5.hit.syb1 <- subset(pks.annotated, SYMBOL %in% syb.1$symbol)
# tn5.hit.syb1.DEG <- subset(DEG, rownames(DEG) %in% tn5.hit.syb1$Geneid)
# tn5.hit.syb1.DEG$group <- "group2"
# 
# 
# tn5.hit.syball.DEG <- rbind(tn5.hit.syb2.DEG, tn5.hit.syb1.DEG)
# 
# 
# ptdata <- reshape2::melt(tn5.hit.syball.DEG[,c(1:3,8)])
# p1 <-
#   ggplot(ptdata, aes(y=value, x=variable, fill=group))+
#   geom_hline(aes(yintercept = 0), color="grey") +
#   geom_boxplot(notch = T, notchwidth = 0.5)+ #
#   themo_demo+labs(title ="log2FC of group1&2 coTGs (Tn5 points)")
#   # ylim(-4,6)
# 
# t.test(tn5.hit.syb2.DEG[,3], tn5.hit.syb1.DEG[,3])
# t.test(tn5.hit.syb2.DEG[,2], tn5.hit.syb1.DEG[,2])
# t.test(tn5.hit.syb2.DEG[,1], tn5.hit.syb1.DEG[,1])
# 
# 
# ########## tn5 counts and coDEGs
# 
# tn5.coTGs.anno <- subset(pks.annotated, pks.annotated$SYMBOL %in% coTGs$target)
# tn5.coTGs.DEG <- subset(DEG, rownames(DEG) %in% tn5.coTGs.anno$Geneid)
# 
# ptdata <- reshape2::melt(tn5.coTGs.DEG[,1:3])
# p2<-
#   ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
#   geom_boxplot(notch = T, notchwidth = 0.5)+
#   themo_demo+labs(title = paste("log2FC distribution of all", 
#                    "NFAT1&RELA co-targtes under 3 stimulations",sep="\n"))+
#   geom_hline(aes(yintercept = 0), color="grey")+
#   ylim(-2,5)
# 
# 
# p3<-
#   ggplot(ptdata, aes(x=value, color=variable)) +
#   geom_density()+
#   xlim(0,4)+ labs(title=paste("log2FC distribution of all",
#                               "NFAT1&RELA co-targtes under 3 stimulations (Tn5 points)",sep="\n"))+
#   themo_demo
# 
# ks.test(tn5.coTGs.DEG[,2], tn5.coTGs.DEG[,3])
# #   
# 
# # # ---------------------------------------------------------------------------
# # # EdgeR
# # # differential expressed gene analysis by edgeR
# # # ---------------------------------------------------------------------------
# 
# library(ggplot2)
# library(reshape2)
# library(ggbeeswarm)
# library("edgeR")
# library(dplyr)
# library(grid)
# library(biomaRt)
# library(org.Hs.eg.db)
# 
# dir_save = "/storage/labnas5/HuangWen5/linwei/"
# setwd(dir_save)
# mycounts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/counts.txt", header = T)
# mycounts2 <- read.table("//storage/labnas5/HuangWen5/linwei/linux/pks_readcounts_bedtools_multicov.txt", header = F)
# colnames(mycounts2) <- c("Chr","Start","End","NC.1","NC.2","I.1","I.2","P.1","P.2","PI.1","PI.2")
# rownames(mycounts2)<- paste("peakID",1:nrow(mycounts2),sep="")
# matrix <- mycounts2[,c(6,7,4,5,8:11)]
# # colnames()
# #分组
# group_list =c(rep("I",2),rep("NC",2),rep("P",2),rep("PI",2))
# 
# counts
# #差异分析
# y <- DGEList(counts=matrix, group = group_list)
# y
# y <- calcNormFactors(y)
# y <- estimateCommonDisp(y)
# y <- estimateTagwiseDisp(y)
# 
# 
# ########################寻找差异基因##############################
# 
# #设计样本matrix
# design <- model.matrix(~0+group, data=y$samples)
# colnames(design) <- levels(y$samples$group)
# 
# fit <- glmQLFit(y, design)
# 
# #一对一比较
# qlf.I_NC  <- glmQLFTest(fit, contrast= c(1,-1,0,0))
# qlf.P_NC <- glmQLFTest(fit, contrast= c(0,-1,1,0))
# qlf.PI_NC <- glmQLFTest(fit, contrast=c(0,-1,0,1))
# 
# #同时比较
# my.contrasts <- makeContrasts(P.NC =P-NC, I.NC =I-NC, PI.NC =PI-NC,
#                               levels=design)
# 
# qlf<- glmQLFTest(fit, contrast=my.contrasts)
# #Example: qlf.S_10minvsNC <- glmQLFTest(fit, contrast=my.contrasts[,"S_10minvsNC"])
# 
# 
# #筛选符合条件的差异基因
# tTag <- topTags(qlf, n=nrow(y))
# # DEG <- subset(tTag$table, PValue < 0.05)
# DEG <- tTag$table
# DEG.pks <- subset(mycounts2, rownames(mycounts2) %in% rownames(DEG))[,1:3]
# DEG.pks$peakID <- rownames(DEG.pks)
# pks <- GRanges(DEG.pks)
# peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
#                              tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
#                              annoDb="org.Hs.eg.db")
# pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)
# 
# 
# pks_readcounts.syb2 <- subset(pks.annotated, SYMBOL %in% syb.2$symbol)
# pks_readcounts.syb2.DEG <- subset(DEG, rownames(DEG) %in% pks_readcounts.syb2$peakID)
# pks_readcounts.syb2.DEG$group <- "group1"
# 
# pks_readcounts.syb1 <- subset(pks.annotated, SYMBOL %in% syb.1$symbol)
# pks_readcounts.syb1.DEG <- subset(DEG, rownames(DEG) %in% pks_readcounts.syb1$peakID)
# pks_readcounts.syb1.DEG$group <- "group2"
# 
# 
# pks_readcounts.syball.DEG <- rbind(pks_readcounts.syb2.DEG, pks_readcounts.syb1.DEG)
# 
# ptdata <- reshape2::melt(pks_readcounts.syball.DEG[,c(1:3,8)])
# p11<-
#   ggplot(ptdata, aes(y=value, x=variable, fill=group))+
#   geom_hline(aes(yintercept = 0), color="grey") +
#   geom_boxplot(notch = T, notchwidth = 0.5)+ #
#   themo_demo+labs(title ="log2FC of group1&2 coTGs")
#   # ylim(-4,6)
# 
# t.test(pks_readcounts.syb2.DEG[,3], pks_readcounts.syb1.DEG[,3])
# t.test(pks_readcounts.syb2.DEG[,2], pks_readcounts.syb1.DEG[,2])
# t.test(pks_readcounts.syb2.DEG[,1], pks_readcounts.syb1.DEG[,1])
# 
# ########## pks counts and coDEGs
# 
# pks.coTGs.anno <- subset(pks.annotated, pks.annotated$SYMBOL %in% coTGs$target)
# pks.coTGs.DEG <- subset(DEG, rownames(DEG) %in% pks.coTGs.anno$peakID)
# 
# ptdata <- reshape2::melt(pks.coTGs.DEG[,1:3])
# p22<-
#   ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
#   geom_boxplot(notch = T, notchwidth = 0.3)+
#   themo_demo+labs(title =paste("log2FC distribution of all",
#                                "NFAT1&RELA co-targtes under 3 stimulations",sep="\n"))+
#   geom_hline(aes(yintercept = 0), color="grey")+
#   ylim(-2,5)
# 
# p33 <-
#   ggplot(ptdata, aes(x=value, color=variable)) +
#   geom_density()+
#   xlim(0,4)+
#   labs(title=paste("log2FC distribution of all",
#                    "NFAT1&RELA co-targtes under 3 stimulations",sep="\n"))+
#   themo_demo
# 
# ks.test(pks.coTGs.DEG[,1], pks.coTGs.DEG[,2])
# ks.test(pks.coTGs.DEG[,2], pks.coTGs.DEG[,3])
# ks.test(pks.coTGs.DEG[,1], pks.coTGs.DEG[,3])
# 
# 
# ### PI up-regulated coTGs
# 
# pks.PIcoTGs.anno <- subset(pks.annotated, pks.annotated$SYMBOL %in% coDEG.PI$Symbol)
# pks.PIcoTGs.DEG <- subset(DEG, rownames(DEG) %in% pks.PIcoTGs.anno$peakID)
# 
# ptdata <- reshape2::melt(pks.PIcoTGs.DEG[,1:3])
# 
# p44 <- 
#   ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
#   geom_boxplot(notch = F, notchwidth = 0.5)+
#   themo_demo+labs(title ="PTS-count-FC of 127 coTGs up-regulated under PMA&Iono")+
#   geom_hline(aes(yintercept = 0), color="grey")
# 
# 
# # # pdf("//storage/labnas5/HuangWen5/linwei/linux/ATAC_seq_2methods.pdf", paper="A4")
# ggpubr::ggarrange(p1,p2,p3, nrow=3, ncol=1)
# ggpubr::ggarrange(p11,p22,p33, nrow=3, ncol=1)
# # dev.off()
# 
# 
# 
# 
# 
# # ############################################ ATAC-seq data ################################################################
# # ##########################        combine heatmap in Fig5 with ATAC
# # ###########################################################################################################################
# # Step1. get promoter region
# ###### ____________________________________________________________________________________________________________________________________________________
library(org.Hs.eg.db)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
columns <- c("tx_name", "gene_id")

PR <- getPromoters(txdb, upstream=5000, downstream=5000, by="gene")
hg38.tss.5kb_bed <- data.frame(PR)
hg38.tss.5kb_bed <- subset(hg38.tss.5kb_bed, hg38.tss.5kb_bed$start > 0 & seqnames %in% c(paste("chr",1:22,sep=""), "chrX","chrY"))
hg38.tss.5kb_bed$peakID <- paste("peakID",1:nrow(hg38.tss.5kb_bed), sep="")
# colnames(hg38.tss.5kb_bed)[1] <- "Chr"
pks <- GRanges(hg38.tss.5kb_bed)
peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                               tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
                               annoDb="org.Hs.eg.db")
pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)

PTS.coTGS <- subset(pks.annotated, pks.annotated$SYMBOL %in% coTGs$target)
PTS.coTGS.bed <- subset(hg38.tss.5kb_bed, hg38.tss.5kb_bed$peakID%in% PTS.coTGS$peakID)

# write.table(PTS.coTGS.bed[,c(1:3,6)], "//storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS.coTGS.bed", sep = "\t",
#                                    row.names = F, col.names = F, quote=F)
# sort -k 1,1 -k2,2n PTS.coTGS.bed > PTS.coTGS_sorted.bed
# //storage/disk2/HuangW/HW_PI_timepoint_ATAC/03.Post_aligned_filtering/last.bam/
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-240-1.last.bam rep2/renamed.bam/I-240-2.last.bam rep1/renamed.bam/P-240-1.last.bam rep2/renamed.bam/P-240-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-240-2.last.bam -bed //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS.coTGS.bed > //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS_coTGs_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-30-1.last.bam rep2/renamed.bam/I-30-2.last.bam rep1/renamed.bam/P-30-1.last.bam rep2/renamed.bam/P-30-2.last.bam rep1/renamed.bam/PI-30-1.last.bam rep2/renamed.bam/PI-30-2.last.bam -bed //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS.coTGS.bed > //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS30min_coTGs_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-120-1.last.bam rep2/renamed.bam/I-120-2.last.bam rep1/renamed.bam/P-120-1.last.bam rep2/renamed.bam/P-120-2.last.bam rep1/renamed.bam/PI-30-1.last.bam rep2/renamed.bam/PI-30-2.last.bam -bed //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS.coTGS.bed > //storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS120min_coTGs_bedtools_multicov.txt
# 
# #################################################### PTS bedtools multicov HW0401 #############################################################
pks_readcounts <- read.table("//storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/PTS/PTS120min_coTGs_bedtools_multicov.txt", header = F)
colnames(pks_readcounts) <- c("chr","start","end","peakID","NC.1","NC.2","I.1","I.2","P.1","P.2","PI.1","PI.2")

X <- data.frame(pks_readcounts[,c(4,1:3)],
                NC= (pks_readcounts$NC.1+pks_readcounts$NC.2)/2,
                I= (pks_readcounts$I.1+pks_readcounts$I.2)/2,
                P= (pks_readcounts$P.1+pks_readcounts$P.2)/2,
                PI= (pks_readcounts$PI.1+pks_readcounts$PI.2)/2)

pks.annotated <- merge(X, PTS.coTGS, by="peakID")

pks.annotated.FC <- data.frame(pks.annotated[,c(1:4,19)],
                                I= pks.annotated$I/(pks.annotated$NC+0.01),
                                P= pks.annotated$P/(pks.annotated$NC+0.01),
                                PI= pks.annotated$PI/(pks.annotated$NC+0.01))
### all coTGs
ptdata <- reshape2::melt(pks.annotated.FC[,c(6:8)],)
b1 <-
  ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  themo_demo+labs(title ="PTS-count-FC of 1345 coTGs")+
  geom_hline(aes(yintercept = 1), color="grey")+ylim(0,2)

### PI up-regulated coTGs
pt <- subset(pks.annotated.FC, SYMBOL%in% coDEG.PI$Symbol)
pt2 <- reshape2::melt(pt[,c(6:8)],)
b2 <-
  ggplot(pt2, aes(y=value, x=variable, fill=variable))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  themo_demo+labs(title ="PTS-count-FC of 115 coTGs up-regulated under PMA&Iono")+
  geom_hline(aes(yintercept = 1), color="grey")+ylim(0,2)


## 2 group coTGs
pks.annotated.FC.syb2 <- subset(pks.annotated.FC, SYMBOL %in% syb.2$symbol)
pks.annotated.FC.syb2$group= "group1"
pks.annotated.FC.syb1 <- subset(pks.annotated.FC, SYMBOL %in% syb.1$symbol)
pks.annotated.FC.syb1$group= "group2"

pks.annotated.FC.syball <- rbind(pks.annotated.FC.syb2, pks.annotated.FC.syb1)

pt2<- reshape2::melt(pks.annotated.FC.syball[,c(6:9)],)
b3 <-
  ggplot(pt2, aes(y=value, x=variable, fill=group))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  themo_demo+labs(title ="PTS-count-FC of 44 coTGs in 2 groups")+
geom_hline(aes(yintercept = 1), color="grey")+
  ylim(0,2)

t.test(subset(pks.annotated.FC.syball, group=="group1")[,9], subset(pks.annotated.FC.syball, group=="group2")[,9])
t.test(subset(pks.annotated.FC.syball, group=="group1")[,8], subset(pks.annotated.FC.syball, group=="group2")[,8])
t.test(subset(pks.annotated.FC.syball, group=="group1")[,7], subset(pks.annotated.FC.syball, group=="group2")[,7])

ggpubr::ggarrange(b1,b2,b3, nrow=1, ncol=3)


########## pks counts and coDEGs
#
# # pks.coTGs.counts <- subset(pks_readcounts, pks_readcounts$peakID %in% pks.coTGs.anno$peakID)
# #
# # pks.coTGs.meancounts <- data.frame(pks.coTGs.counts[,c(1:3,12)],
# #                                    NC=(pks.coTGs.counts$V4+pks.coTGs.counts$V5)/2,
# #                                    I=(pks.coTGs.counts$V6+pks.coTGs.counts$V7)/2,
# #                                    P=(pks.coTGs.counts$V8+pks.coTGs.counts$V9)/2,
# #                                    PI=(pks.coTGs.counts$V10+pks.coTGs.counts$V11)/2)
# #
# # pks.coTGs.FC <- data.frame(tn5.coTGs.meancounts[,1:4],
# #                            I= tn5.coTGs.meancounts$I/(tn5.coTGs.meancounts$NC+0.001),
# #                            P= tn5.coTGs.meancounts$P/(tn5.coTGs.meancounts$NC+0.001),
# #                            PI= tn5.coTGs.meancounts$PI/(tn5.coTGs.meancounts$NC+0.001))
# #
# # ptdata <- reshape2::melt(pks.coTGs.FC[,5:7])
# # ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
# #   geom_boxplot(notch = T, notchwidth = 0.5)+
# #   themo_demo+labs(title ="Fold change of peak counts")+
# #   geom_hline(aes(yintercept = 1), color="grey")+
# #   ylim(0,5)
# #
# #   ggplot(ptdata, aes(x=value, color=variable)) +
# #   geom_density()+
# #   xlim(0,4)+ labs(title=paste("Peak-count-FC distribution of all",
# #                               "NFAT1&RELA co-targtes under 3 stimulations",sep="\n"))+
# #   themo_demo
# #
# # ks.test(pks.coTGs.FC[,5], pks.coTGs.FC[,6])
# # ks.test(pks.coTGs.FC[,5], pks.coTGs.FC[,7])
# # ks.test(pks.coTGs.FC[,6], pks.coTGs.FC[,7])
#
# # #################################################### Back2peaks bedtools multicov HW0401 #############################################################
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-30-1.last.bam rep2/renamed.bam/I-30-2.last.bam rep1/renamed.bam/P-30-1.last.bam rep2/renamed.bam/P-30-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-30-2.last.bam -bed //storage/labnas5/HuangWen5/linwei/linux/pks_in_counts_hw.bed > //storage/labnas5/HuangWen5/linwei/pks30min_readcounts_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-120-1.last.bam rep2/renamed.bam/I-120-2.last.bam rep1/renamed.bam/P-120-1.last.bam rep2/renamed.bam/P-120-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-120-2.last.bam -bed //storage/labnas5/HuangWen5/linwei/linux/pks_in_counts_hw.bed > //storage/labnas5/HuangWen5/linwei/pks120min_readcounts_bedtools_multicov.txt
# bedtools multicov -bams rep1/renamed.bam/NC-0-1.last.bam rep2/renamed.bam/NC-0-2.last.bam rep1/renamed.bam/I-240-1.last.bam rep2/renamed.bam/I-240-2.last.bam rep1/renamed.bam/P-240-1.last.bam rep2/renamed.bam/P-240-2.last.bam rep1/renamed.bam/PI-240-1.last.bam rep2/renamed.bam/PI-240-2.last.bam -bed //storage/labnas5/HuangWen5/linwei/linux/pks_in_counts_hw.bed > //storage/labnas5/HuangWen5/linwei/pks240min_readcounts_bedtools_multicov.txt
# 

# pks_readcounts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/pks30min_readcounts_bedtools_multicov.txt", header = F)
pks_readcounts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/pks120min_readcounts_bedtools_multicov.txt", header = F)
# pks_readcounts <- read.table("//storage/labnas5/HuangWen5/linwei/linux/pks240min_readcounts_bedtools_multicov.txt", header = F)


colnames(pks_readcounts) <- c("chr","start","end","NC.1","NC.2","I.1","I.2","P.1","P.2","PI.1","PI.2")
pks_readcounts$peakID <- paste("peak",1:nrow(pks_readcounts), paste="")

X <- data.frame(pks_readcounts[,c(12,1:3)],
                NC= (pks_readcounts$NC.1+pks_readcounts$NC.2)/2,
                I= (pks_readcounts$I.1+pks_readcounts$I.2)/2,
                P= (pks_readcounts$P.1+pks_readcounts$P.2)/2,
                PI= (pks_readcounts$PI.1+pks_readcounts$PI.2)/2)


pks <- GRanges(X[,1:4])
peakAnnoList <- annotatePeak(pks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                             tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,
                             annoDb="org.Hs.eg.db")
pks.annotated <- data.frame(peakAnnoList@anno@elementMetadata@listData)
pks.annotated <- merge(X, pks.annotated, by="peakID")

pks.annotated.FC <- data.frame(pks.annotated[,1:4], I= pks.annotated$I/(pks.annotated$NC+0.1), 
                               P= pks.annotated$I/(pks.annotated$NC+0.1),  
                               PI= pks.annotated$PI/(pks.annotated$NC+0.1), pks.annotated[,c(9,17,19,20)] )
symbol.for.pks <- data.frame(symbol= unique(pks.annotated$SYMBOL))

max.pks.annotated.FC <- data.frame()
for (syb in 1:nrow(symbol.for.pks)){
 x <- subset(pks.annotated.FC, SYMBOL == symbol.for.pks[syb,1])
 max.pks.annotated.FC <- rbind(max.pks.annotated.FC, x[which.max(x$PI),])
}


### all coTGs
ptdata <- subset(max.pks.annotated.FC, SYMBOL%in% coTGs$target)
ptdata <- reshape2::melt(ptdata[,c(5:7)],)
b1 <-
  ggplot(ptdata, aes(y=value, x=variable, fill=variable))+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  themo_demo+labs(title ="PTS-count-FC of 1376 coTGs")+
  geom_hline(aes(yintercept = 1), color="grey")+ylim(0,4)

### PI up-regulated coTGs
pt <- subset(max.pks.annotated.FC, SYMBOL%in% coDEG.PI$Symbol)
pt2 <- reshape2::melt(pt[,c(5:7)],)
b2 <-
  ggplot(pt2, aes(y=as.numeric(value), x=variable, fill=variable))+
  geom_signif(comparisons =COMP ,step_increase = 0.1, map_signif_level = F,test = t.test)+
  geom_boxplot(notch = F, notchwidth = 0.5)+
  themo_demo+labs(title ="PTS-count-FC of 127 coTGs up-regulated under PMA&Iono")+
  geom_hline(aes(yintercept = 1), color="grey")+ylim(0,4)


## 2 group coTGs
pks.annotated.FC.syb2 <- subset(max.pks.annotated.FC, SYMBOL %in% syb.2$symbol)
pks.annotated.FC.syb2$group= "group1"
pks.annotated.FC.syb1 <- subset(max.pks.annotated.FC, SYMBOL %in% syb.1$symbol)
pks.annotated.FC.syb1$group= "group2"

pks.annotated.FC.syball <- rbind(pks.annotated.FC.syb2, pks.annotated.FC.syb1)

pt2<- reshape2::melt(pks.annotated.FC.syball[,c(5:7,12)],)
b3 <-
  ggplot(pt2, aes(y=value, x=variable, fill=group))+
  geom_boxplot(notch = F, notchwidth = 0.5)+  
  geom_signif(comparisons =COMP ,step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo+labs(title ="PTS-count-FC of 51 coTGs in 2 groups")+
  geom_hline(aes(yintercept = 1), color="grey")+
  ylim(0,4)

t.test(subset(pks.annotated.FC.syball, group=="group1")[,7], subset(pks.annotated.FC.syball, group=="group2")[,7])
t.test(subset(pks.annotated.FC.syball, group=="group1")[,6], subset(pks.annotated.FC.syball, group=="group2")[,6])
t.test(subset(pks.annotated.FC.syball, group=="group1")[,5], subset(pks.annotated.FC.syball, group=="group2")[,5])

# g1 <- ggpubr::ggarrange(b1,b2,b3, nrow=1, ncol=3)
g2 <- ggpubr::ggarrange(b1,b2,b3, nrow=1, ncol=3)
# g3 <- ggpubr::ggarrange(b1,b2,b3, nrow=1, ncol=3)

ggpubr::ggarrange(g1,g2,g3, nrow=3, ncol=1)

