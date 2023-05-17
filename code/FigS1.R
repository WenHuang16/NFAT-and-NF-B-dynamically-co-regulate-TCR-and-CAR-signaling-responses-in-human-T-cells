################################################ Step1 ###############################################################################
library(dplyr)

in.dir <- "//storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/"
setwd(in.dir)

TF_TG <- dir(paste(in.dir, "/5kb/", sep=""))
TF_TG.index <- data.frame(TF_TG)

celltypeList <- read.table(paste(in.dir,"celltypeList.tab", sep=""), sep="\t",header = T)

#### get all experiments ID conducted in CD4+ T
CD4T.ID.in.celltypeList <- subset(celltypeList,
                                  Cell_type %in% c("Th1 Cells", "CD4+ Th1","CD4-Positive T-Lymphocytes", "CD4+") &
                                  Genome %in% c("hg38"))

CD4TID.all <- c()
for (r in 1:nrow(CD4T.ID.in.celltypeList)){
  IDall <- CD4T.ID.in.celltypeList$ID[r]
  CD4TID <- strsplit(IDall, "[,]")
  CD4TID.all <- c(CD4TID.all,c(CD4TID[[1]]))
}


#### get new TF-TG table contain experiments conducted in blood cell-lines
TF_TG_CD4T <- list()
kkk=0
for (tfi in 1:length(TF_TG)){
  ### tfi = 546

  TF.table <- read.csv(paste(in.dir, "/5kb/",TF_TG[tfi], sep=""), sep="\t")
  ID <- colnames(TF.table)[3:(ncol(TF.table)-1)]

  y <- strsplit(ID, "[.]")

  TF.ID <- data.frame(ID=NA, cellline=NA)
  selected.ID <- data.frame()
  for (i in 1:length(y)){
    TF.ID <- na.omit(rbind(TF.ID, c(ID = y[[i]][1], cellline = y[[i]][2])) )
    if (TF.ID$ID[i] %in% CD4TID.all){
      selected.ID[i,1] <- 1
    }else{
      selected.ID[i,1]  <- 0
    }
  }

  TF.CD4T.id <- as.numeric(rownames(subset(selected.ID, selected.ID$V1 == 1)))+2
  TF.table.CD4T <- TF.table[,c(1,2,TF.CD4T.id)]

  if (ncol(TF.table.CD4T) > 2){ #TF without CD4T samples will be filtered due to a ncol(TF.table.CD4T)=2
    kkk = kkk+1
    TF.table.CD4T$CD4T.Average <- apply(TF.table.CD4T[3:ncol(TF.table.CD4T)], 1, mean )
    # threshold <- data.frame(quantile(TF.table.CD4T$CD4T.Average)) ### threshold was not used
    TF.table.CD4T<- subset(TF.table.CD4T, TF.table.CD4T$CD4T.Average > 0)#threshold[1,1]) #### 0,20

    TF_TG_CD4T[[kkk]] <- TF.table.CD4T
    namei <- strsplit(TF_TG[tfi], "[.]")
    names(TF_TG_CD4T)[kkk] <- namei[[1]][1]
  }

}

############################################ Get putative target genes 3hubs #######################################################################

MAPK <- c("JUN","JUN","JUND","FOS","FOSB","FOSL1","FOSL2")
Ca <- c("NFATC1","NFATC2","NFATC3")
PKC <- c("RELA","REL", "RELB", "NFKB1","NFKB2")

hubs <- c("MAPK","Ca","PKC")
Len = "5kb"#c("1kb", "5kb" ,"10kb")


## for (ri in c(1,5,10)){
ri=5
Li <- paste(ri,"kb", sep="")

for (i in 1:length(hubs)){
  hub <- get(hubs[i])
  TGs <- data.frame(matrix(0,1,2))
  colnames(TGs) <- c("target", "CD4T.norm")
  for (k in 1:length(hub)){
    # TGik <- read.csv(paste(in.dir, "/",Li,"/",hub[k],".",ri,".tsv", sep=""), sep="\t")
    # TGik <- subset(TGik, TGik[,2] > 20)

    TGik <- TF_TG_CD4T[[hub[k]]]
    TGik$CD4T.norm <- (TGik$CD4T.Average-min(TGik$CD4T.Average))/(max(TGik$CD4T.Average)-min(TGik$CD4T.Average))# max-min normalization(score)
    #TGs <- unique(c(TGs,TGik$Target_genes))
    TGs <- rbind(TGs, data.frame(target= TGik$Target_genes, CD4T.norm= TGik$CD4T.norm))
  }

  TGs.unique <- aggregate(CD4T.norm ~ target, data = TGs, sum)
  TGs.unique <- subset(TGs.unique, !TGs.unique$target == "0")
  assign(paste("TG.",hubs[i],sep=""),TGs.unique)
}

coTGs <- intersect(TG.Ca$target, TG.PKC$target )
coTGs.1 <- subset(TG.Ca, TG.Ca$target %in% coTGs)
coTGs.2 <- subset(TG.PKC, TG.PKC$target %in% coTGs)
coTGs <- merge(coTGs.1, coTGs.2, by="target")
colnames(coTGs) <- c("target","NFAT.score","NFKB.score")
coTGs$co.score<- coTGs$NFAT.score + coTGs$NFKB.score

coTGs.hi <- subset(coTGs, coTGs$co.score > as.numeric(quantile(coTGs$co.score)[4]))
# 
# 
# write.csv(coTGs, "/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/coTGs.all_2022_Home.csv", row.names = F )
# write.csv(coTGs.hi,"/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/coTGs.top25_2022_Home.csv", row.names = F )
# 
setwd("/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/")

library(VennDiagram)
# venn.plot <- venn.diagram(x = list(MAPK = TG.MAPK$target, Ca = TG.Ca$target, PKC = TG.PKC$target),
#                           filename = paste("Venn.3hubs_",Li, "_2022_Home.tiff",sep=""), imagetype = "tiff",
#                           lwd = 3, fill = c("cornflowerblue", "darkorchid1", "darkgreen"),
#                           alpha = 0.6, label.col = "white", cex = 1.5,
#                           fontfamily = "serif", fontface = "bold",
#                           # cat.col = c("cornflowerblue", "darkorchid1"),
#                           cat.cex = 2,
#                           cat.fontfamily = "serif",
#                           cat.fontface = "bold",
#                           margin = 0.05,
#                           cat.dist = c(0.03, 0.03, 0.03),
#                           # cat.pos = c(-20, 20)
# )
# y <- calculate.overlap(x = list(MAPK = TG.MAPK$target, Ca = TG.Ca$target, PKC = TG.PKC$target))
# pure.MAPK <- data.frame(symbol= y[["a1"]])
# pure.PKC <- data.frame(symbol=y[["a7"]])
# pure.Ca <- data.frame(symbol=y[["a3"]])
#
# overlap.3 <- data.frame(symbol=y[["a5"]])
# overlap.MAPK_PKC <- data.frame(symbol=y[["a4"]])
# overlap.MAPK_Ca <- data.frame(symbol=y[["a2"]])
# overlap.PKC_Ca <- unique(data.frame(symbol=c(y[["a5"]], y[["a6"]])))
#
# coTGs <- data.frame(symbol = c(overlap.3$symbol, overlap.PKC_Ca$symbol))
#
# write.csv(pure.MAPK, paste("pure.MAPK.",Li,"_2022_Home.csv",sep=""), row.names = F)
# write.csv(pure.PKC, paste("pure.PKC.",Li,"_2022_Home.csv",sep=""), row.names = F)
# write.csv(pure.Ca, paste("pure.Ca.",Li,"_2022_Home.csv",sep=""), row.names = F)
# write.csv(overlap.3, paste("overlap.3.",Li,"_2022_Home.csv",sep=""), row.names = F)
# write.csv(overlap.PKC_Ca, paste("overlap.PKC_Ca",Li,"_2022_Home.csv",sep=""), row.names = F)
#
#
#


# save.image("/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/CD4_TF_TG_analysis_Step1_2023.RData")

################################################# Step1 ###############################################################################
load("/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/CD4_TF_TG_analysis_Step1_2023.RData")




################################################# Step2 ###############################################################################


################# NFAT1-RELA CO-REGULATED GENE
################# Add CO-REGULATED GENE inferred by dorothea
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
TG.NFAT.dorothea <- subset(dorothea_hs, dorothea_hs$tf %in% c("NFATC1","NFATC2","NFATC3") & dorothea_hs$confidence %in% c("A","B","C"))
TG.NFKB.dorothea <- subset(dorothea_hs, dorothea_hs$tf %in% c("RELA","RELB","REL","NFKB1","NFKB2") & dorothea_hs$confidence %in% c("A","B","C"))
coDoR.TG <- intersect(TG.NFAT.dorothea$target, TG.NFKB.dorothea$target)
# write.csv(unique(c(coTGs$target,coDoR.TG)), "/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/coTGs.all.addDorothea_2022_Home.csv", row.names = F )



################# NFAT1-RELA CO-REGULATED GENE

GO.typ <- readxl::read_xlsx("/storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/coTGs_all.tqy_j6yub/metascape_result.xlsx", sheet = "Enrichment")

GO.typ <- GO.typ[grep("Summary",GO.typ$GroupID),]
GO.typ <- GO.typ[grep("GO Bio",GO.typ$Category),]
GO.typ$Genenumber <- NA

for (nr in 1:nrow(GO.typ)){
  GO.typ$Genenumber[nr] = length(strsplit(GO.typ$Genes[nr], split=",")[[1]])
}

nr <- 10
go_enrich_df<-data.frame(ID=GO.typ$GroupID[1:nr],
                         Description=GO.typ$Description[1:nr],
                         `-logP`= -(GO.typ$LogP)[1:nr],
                         type=GO.typ$Category[1:nr],
                         Genenumber = GO.typ$Genenumber[1:nr])

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## add GO terms description names
labels <- go_enrich_df$Description
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar
CPCOLS <- "#eb6bac"

p.GO <- ggplot(data=go_enrich_df, aes(x=number, y=X.logP, fill=Genenumber)) +coord_flip() +
  geom_bar(stat="identity", width=0.5) +
  # scale_fill_manual(values = CPCOLS) +
  scale_fill_gradient(low="LightSkyBlue2",high="orangeRed")+
  theme_test() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = paste("GO analysis of NFAT/NFkB 1548Entrez co-targets", sep="")) #barplot_coTGs.1548_Entrez ##barplot_coTGs.1557_Entrez

# ggsave("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/coTGs.all_PKC_Ca_GO_barplot.pdf", width = 10, height = 8)

################# CO-REGULATED GENE "WEIGHT" of TF PAIRS 

TF_weights <- data.frame()
system.time(
  for (i in 1:(length(TF_TG_CD4T)-1)){
    for (j in (i+1):length(TF_TG_CD4T)){
      x <- data.frame(TF1 = names(TF_TG_CD4T)[i],
                      TF2 = names(TF_TG_CD4T)[j],
                      weight = length(intersect(TF_TG_CD4T[[i]]$Target_genes, TF_TG_CD4T[[j]]$Target_genes)))
      
      TF_weights <- rbind(TF_weights, x)
      print(c(i,j))
    }
  })
# 
## user  system elapsed
## 137.471   2.300 138.666
# 
TF_min_ratio <- data.frame()

for (i in 1:(length(TF_TG_CD4T)-1)){
  for (j in (i+1):length(TF_TG_CD4T)){
    y <- data.frame(TF1.num =nrow(TF_TG_CD4T[[i]]), TF2.num = nrow(TF_TG_CD4T[[j]]))
    TF_min_ratio <- rbind(TF_min_ratio, y)
    
    print(c(i,j))
  }
}


TF_weights$TF1.ratio <- TF_weights$weight/TF_min_ratio$TF1.num
TF_weights$TF2.ratio <- TF_weights$weight/TF_min_ratio$TF2.num
TF_weights$pair <- paste(TF_weights$TF1,TF_weights$TF2,sep=".")



############################################ hypergeometric test ############################################

gn <- 25000 ######
# 

kk <- cbind(TF_weights, TF_min_ratio)
kk$prob.random <- (kk$TF1.num/gn) * (kk$TF2.num/gn)
kk$prob.real <- kk$weight/gn
kk$prob.FC <- kk$prob.real/kk$prob.random

for (xx in 1:nrow(kk)){
  weight = kk$weight[xx]
  TF1.num = kk$TF1.num[xx]
  TF2.num = kk$TF2.num[xx]
  
  kk$p.value[xx] <- phyper(weight-1, TF1.num, gn-TF1.num, TF2.num, lower.tail=F, log=F)
  kk$`-logP`[xx] <- -phyper(weight-1, TF1.num, gn-TF1.num, TF2.num, lower.tail=F, log=T) ## ln
  
  print(xx)
}

# kk$p.adj <- p.adjust(exp(kk$`-logP`), method = "BH", n = nrow(kk))
kk$p.adj <- p.adjust(kk$p.value, method = "BH", n = nrow(kk))


# ### p-value matrix #####
# ##
pvalue.TF_TG_CD4T <- data.frame(matrix(1, nrow=length(TF_TG_CD4T),ncol= length(TF_TG_CD4T)))
colnames(pvalue.TF_TG_CD4T) = rownames(pvalue.TF_TG_CD4T) = names(TF_TG_CD4T)
K=length(TF_TG_CD4T)

for (i in 1:(K-1)){
  for (j in (i+1):K){
    p <- subset(kk, kk$TF1 == colnames(pvalue.TF_TG_CD4T)[i] & kk$TF2 == rownames(pvalue.TF_TG_CD4T)[j])
    pvalue.TF_TG_CD4T[i,j] = pvalue.TF_TG_CD4T[j,i] = as.numeric(p$`-logP`)
  }
}

# dist(pvalue.TF_TG_CD4T)

# pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1B_TF_TG_-logP_heatmap_CD4T_2023.pdf", width = 10, height = 10)
library(pheatmap)
# display.brewer.all()
pheatmap(pvalue.TF_TG_CD4T, clustering_method = "ward.D",
         #clustering_distance_cols = "euclidean",
         treeheight_row = 0, treeheight_col = 60,
         # color = colorRampPalette(brewer.pal(9,"Blues"))(100),
         color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
         # legend_breaks = 0:16,
         # legend_labels = c("0","8","16"),
         cutree_cols = 3,cutree_rows = 3,
         fontsize = 10)

# dev.off()


library("ggplot2")
library("RColorBrewer")
library("ggrepel")
#############

den.plot <- rbind(data.frame(ID = rep("All", nrow(kk)), weight = kk$weight ),
                  data.frame(ID = rep("pPass", nrow(plot.TF_weights_pPass)), weight = plot.TF_weights_pPass$weight ))
density(log(kk$weight+1))
# 
# ggplot(kk,  aes(x= log(weight+1))) +
#   
#   geom_vline(xintercept = log(1531+1), colour = "red")+ #NFATC2.RELA
#   geom_vline(xintercept = log(1085+1), colour = "pink")+ # Median
#   geom_vline(xintercept = quantile(log(kk$weight+1))[[2]] )+ # Median
#   geom_vline(xintercept = quantile(log(kk$weight+1))[[3]] )+ # Median
#   geom_vline(xintercept = quantile(log(kk$weight+1))[[4]] )+ # 3rd Qu.
#   geom_vline(xintercept =  log(694+1), colour = "orange")+ # FOS.JUN
#   geom_vline(xintercept =  log(142+1), colour = "red", linetype="dashed")+ # NFATC1.RELA
#   labs(title = "red:NFATC2.RELA,orange:FOS.JUN")+
#   geom_density(color = "black", fill = "gray")+
#   xlim(-1,11)+
#   ylim(0,0.2)+
#   theme_bw()+
#   theme(panel.grid.major =element_blank(),  # 鍘婚櫎杈规
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_text(color='black'),
#         plot.title = element_text(hjust = 0.5),
#         title=element_text(size = 15),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=20))
# 
# dev.off()
kkord <- kk[order(kk$weight, decreasing = T),]
kkord$rank <- 1:nrow(kkord)
84/703

############## p.adj < 0.05, Observed-to-random ratio #########################################################
# pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1C_TF_TG_logWeight_density_CD4T_2023.pdf", width = 10, height = 6)

kk2 <- subset(kk, kk$p.adj < 0.05)

ggplot(kk2,  aes(x= prob.FC)) +
  geom_density(color = "black", fill = "gray")+
  geom_vline(xintercept = 2.240517, colour = "red")+ #NFATC2.RELA
  geom_vline(xintercept = 2.39, colour = "red")+ #NFKB2.RELA
  geom_vline(xintercept = quantile(kk2$prob.FC)[[2]] )+ # 1rd 
  geom_vline(xintercept = quantile(kk2$prob.FC)[[3]] )+ # Median
  geom_vline(xintercept = quantile(kk2$prob.FC)[[4]] )+ # 3rd Qu.
  geom_vline(xintercept =  6.142497, colour = "orange")+ # FOS.JUN
  geom_vline(xintercept =  2.301089, colour = "red", linetype="dashed")+ # NFATC1.RELA
  xlim(0,10)+
  labs(title = "red:NFATC2.RELA,orange:FOS.JUN, p.adj<0.05")+
  theme_bw()+
  theme(panel.grid.major =element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(color='black'),
        plot.title = element_text(hjust = 0.5),
        title=element_text(size = 15),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))

dev.off()

############## Ranking of all possible TF pairs by weight #########################################################
pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1G_TF_TG_logWeight_density_CD4T_2023.pdf", width = 10, height = 6)

# den.plot <- rbind(data.frame(ID = rep("All", nrow(kk)), weight = kk$weight ),
#                   data.frame(ID = rep("pPass", nrow(plot.TF_weights_pPass)), weight = plot.TF_weights_pPass$weight ))
density(log(kk$weight+1))

ggplot(kk,  aes(x= log2(weight+1))) +
  geom_density(color = "black", fill = "gray")+
  geom_vline(xintercept = log(1531+1), colour = "red")+ #NFATC2.RELA
  geom_vline(xintercept = quantile(log(kk$weight+1))[[3]] )+ # Median
  geom_vline(xintercept = quantile(log(kk$weight+1))[[4]] )+ # 3rd Qu.
  geom_vline(xintercept =  log(694+1), colour = "orange")+ # FOS.JUN
  geom_vline(xintercept =  log(142+1), colour = "red", linetype="dashed")+ # NFATC1.RELA
  labs(title = "red:NFATC2.RELA,orange:FOS.JUN")+
  geom_density(color = "black", fill = "gray")+
  xlim(0,15)+
  ylim(0,0.2)+
  theme_bw()+
  theme(panel.grid.major =element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(color='black'),
        plot.title = element_text(hjust = 0.5),
        title=element_text(size = 25),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))

dev.off()

############################################# co-number plot ############################################
plot.TF_weights <- kk[order(-kk$weight), ]
plot.TF_weights$pair <- factor(plot.TF_weights$pair, levels=plot.TF_weights$pair,ordered = TRUE)
plot.TF_weights$NO. <- c(1:nrow(plot.TF_weights))

p<- ggplot(plot.TF_weights,  aes(x=`NO.`, y= log(weight+1)))+
  geom_line() +
  # geom_area(fill = "SkyBlue", alpha = 0.5)+
  theme_bw()+
  theme(panel.grid.major =element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(color='black'),
        plot.title = element_text(hjust = 0.5),
        title=element_text(size = 25),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))


# subset(plot.TF_weights, plot.TF_weights$TF1.ratio <0.8 & plot.TF_weights$TF2.ratio <0.8 & (plot.TF_weights$TF1 == "RELA" |plot.TF_weights$TF1 == "RELA" )) 

labelx <- subset(plot.TF_weights, plot.TF_weights$pair %in% c("FOS.JUNB","NFATC2.RELA", "NFATC1.RELA","NFKB1.RELA"))

p1 <- p + geom_text_repel(data=labelx,  aes(x=`NO.` , y= log(weight+1), label= pair),
                          lineheight=3, size=5, color="Blue",alpha=0.9,fontface='italic',
                          box.padding=unit(0.5,"lines"),
                          segment.color='red',segment.alpha=0.5,segment.size=1)  

p2 <- p1 + #geom_hline(yintercept = log2(1707+1)) + geom_vline(xintercept = 28034) +
  # geom_vline(xintercept = nrow(plot.TF_weights)*0.5,  color="grey")+
  geom_hline(yintercept = as.numeric(quantile(log(plot.TF_weights$weight+1))[2]), color="grey")+
  geom_hline(yintercept = as.numeric(quantile(log(plot.TF_weights$weight+1))[3]), color="grey")+
  geom_hline(yintercept = as.numeric(quantile(log(plot.TF_weights$weight+1))[4]), color="grey")
  
# pdf("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1G_TF_TG_rank_logWeight_density_CD4T_2023.pdf", width = 10, height = 5)
p2
p
# dev.off()

# svg(file="//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1G_TF_TG_rank_logWeight_density_CD4T_2023.svg", width = 10, height = 5)
p2
p
# dev.off()

# ################## subtract co-regulated genes for GO analysis 0916 ##3
# TG_NFATC1 <- TF_TG_CD4T[["NFATC1"]]
TG_NFATC2 <- TF_TG_CD4T[["NFATC2"]]
TG_RELA <- TF_TG_CD4T[["RELA"]]
CoTGs.NFAT1_RELA <- unique(intersect(TG_NFATC2$Target_genes, TG_RELA$Target_genes))
# CoTGs.NFAT2_RELA <- unique(intersect(TG_NFATC1$Target_genes, TG_RELA$Target_genes))
# 
##### Not co-regulated genes!!!
onlyTGs.NFAT1 <- unique(subset(TG_NFATC2, !Target_genes %in% TG_RELA$Target_genes)$Target_genes)
onlyTGs.RELA <- unique(subset(TG_RELA, !Target_genes %in% TG_NFATC2$Target_genes)$Target_genes)

# write.csv(onlyTGs.NFAT1, "//storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/onlyTGs.NFAT1.csv", row.names = F )
# write.csv(onlyTGs.RELA[1:1000],"//storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/onlyTGs.RELA_top1000.csv", row.names = F )
# write.csv(onlyTGs.RELA[1:300],"//storage/labnas5/HuangWen2023/T cell dynamics paper/TFs_chip-atlas/metadata/onlyTGs.RELA_top300.csv", row.names = F )

# ######### onlyTGs.RELA_top1000 GENE

# GO.typ <- readxl::read_xlsx("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/onlyTGs.RELA_top1000_all.t37eq2erw/metascape_result.xlsx", sheet = "Enrichment")
GO.typ <- readxl::read_xlsx("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/onlyTGs.RELA_top300_all.t_xdyclzs/metascape_result.xlsx", sheet = "Enrichment")
# GO.typ <- readxl::read_xlsx("/storage/disk2/HuangW/T_Cell_Dynamic_paper/metadata/onlyTGs.NFAT1_all.t7b31wmpu/metascape_result.xlsx", sheet = "Enrichment")

GO.typ <- GO.typ[grep("Summary",GO.typ$GroupID),]
GO.typ <- GO.typ[grep("GO Bio",GO.typ$Category),]
GO.typ$Genenumber <- NA
for (nr in 1:nrow(GO.typ)){
  GO.typ$Genenumber[nr] = length(strsplit(GO.typ$Genes[nr], split=",")[[1]])
}

nr <- 10
go_enrich_df<-data.frame(ID=GO.typ$GroupID[1:nr],
                         Description=GO.typ$Description[1:nr],
                         `-logP`= -(GO.typ$LogP)[1:nr],
                         type=GO.typ$Category[1:nr],
                         Genenumber = GO.typ$Genenumber[1:nr])

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## add GO terms description names
labels <- go_enrich_df$Description
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar
CPCOLS <- "#eb6bac"

p.GO <- ggplot(data=go_enrich_df, aes(x=number, y=X.logP, fill=Genenumber)) +coord_flip() +
  geom_bar(stat="identity", width=0.5) +
  # scale_fill_manual(values = CPCOLS) +
  scale_fill_gradient(low="LightSkyBlue2",high="orangeRed")+
  theme_test() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = paste("onlyTGs.RELA_top300", sep=""))
  # labs(title = paste("onlyTGs.RELA_top1000", sep=""))
# 
# ggsave("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1E_1_onlyTGs.NFAT1_barplot_10terms.pdf", width = 10, height = 8)
# ggsave("//storage/labnas5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS1/FigS1E_2_onlyTGs.RELA_top300_barplot_10terms.pdf", width = 10, height = 8)
