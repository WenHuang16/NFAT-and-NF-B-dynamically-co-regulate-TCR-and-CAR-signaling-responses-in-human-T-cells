##__________________________________________________________________________________________________________________
## SINGLE-CELL DATA WERE USED AFTER DROPPING OUT OUTLIERS & SMOOTHING
##__________________________________________________________________________________________________________________

library("R.matlab")
library("ggplot2")

themo_demo= theme_bw()+
  theme(
    text = element_text(size=10),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text=element_text(color='black'),
    plot.title = element_text(hjust = 0.5),
    title=element_text(size = 15),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10))

###### PMA & Iono ________________________________________________________________________________________________________
###### combine 202303 and 202301 data

### 202303
input_dir <- "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/03.TCR+PMATCR+Iono/metadata/202303_scdata/smoothed_NFAT_NFKB/"
allmat <- dir(input_dir)

cellid_all <- data.frame()
NFAT_smooth_all <- data.frame()
NFKB_smooth_all <- data.frame()
for (i in 1:2){
  cellid <- readMat(paste(input_dir, allmat[i], sep=""))
  NFAT_EXCEL_smooth <- readMat(paste(input_dir, allmat[i+2], sep=""))
  NFKB_EXCEL_smooth <- readMat(paste(input_dir, allmat[i+4], sep=""))
  
  cellid_all <- rbind(cellid_all, t(data.frame(cellid)))
  NFAT_smooth_all <- rbind(NFAT_smooth_all, NFAT_EXCEL_smooth[["NFAT.EXCEL.smooth"]])
  NFKB_smooth_all <- rbind(NFKB_smooth_all, NFKB_EXCEL_smooth[["NFKB.EXCEL.smooth"]])
  
}

treatment_id_all <- data.frame(matrix(0,nrow(cellid_all),3))
for (j in 1:nrow(cellid_all)){
  namei <- strsplit(cellid_all$folder[[j]], "\\", fixed=TRUE)
  treatment_id_all[j,2:3 ] <- cellid_all[j,1:2]
  
  if ( namei[[1]][8] %in% c(paste("pos",0:3, sep=""),paste("pos",21:23, sep=""))){
    treatment_id_all[j, 1] <- "Ab" 
  }else if( namei[[1]][8] %in% c(paste("pos",4:7, sep=""))){
    treatment_id_all[j, 1] <- "I0.6+Ab" 
  }else if( namei[[1]][8] %in% c(paste("pos",8:11, sep=""))){
    treatment_id_all[j, 1] <- "P1+Ab" 
  }else if( namei[[1]][8] %in% c(paste("pos",12:15, sep=""))){
    treatment_id_all[j, 1] <- "I0.6" 
  }else if( namei[[1]][8] %in% c(paste("pos",16:19, sep=""))){
    treatment_id_all[j, 1] <- "P1" }
}
colnames(treatment_id_all) <- c("treatment", "cell", "location")


NFAT_with_ID_part1 <- cbind(NFAT_smooth_all, treatment_id_all$treatment)
NFKB_with_ID_part1 <- cbind(NFKB_smooth_all, treatment_id_all$treatment)


#### 202301

input_dir <- "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/03.TCR+PMATCR+Iono/metadata/202301_scdata/smoothed_NFAT_NFKB/"
allmat <- dir(input_dir)

cellid_all <- data.frame()
NFAT_smooth_all <- data.frame()
NFKB_smooth_all <- data.frame()
for (i in 1:3){
  cellid <- readMat(paste(input_dir, allmat[i], sep=""))
  NFAT_EXCEL_smooth <- readMat(paste(input_dir, allmat[i+3], sep=""))
  NFKB_EXCEL_smooth <- readMat(paste(input_dir, allmat[i+6], sep=""))
  
  cellid_all <- rbind(cellid_all, t(data.frame(cellid)))
  NFAT_smooth_all <- rbind(NFAT_smooth_all, NFAT_EXCEL_smooth[["NFAT.EXCEL.smooth"]])
  NFKB_smooth_all <- rbind(NFKB_smooth_all, NFKB_EXCEL_smooth[["NFKB.EXCEL.smooth"]])
  
}

treatment_id_all <- data.frame(matrix(0,nrow(cellid_all),3))
for (j in 1:nrow(cellid_all)){
  namei <- strsplit(cellid_all$folder[[j]], "\\", fixed=TRUE)
  treatment_id_all[j,2:3 ] <- cellid_all[j,1:2]
  
  if (namei[[1]][7] == "2023-01-09" & namei[[1]][9] %in% c(paste("pos",0:4, sep=""))){
    treatment_id_all[j, 1] <- "P1+Ab"
    }else if(namei[[1]][7] == "2023-01-09" & namei[[1]][9] %in% c(paste("pos",15:19, sep=""))){
    treatment_id_all[j, 1] <- "P1"}
}
colnames(treatment_id_all) <- c("treatment", "cell", "location")

#####
NFAT_with_ID_part2 <- cbind(NFAT_smooth_all, treatment_id_all$treatment)
NFKB_with_ID_part2 <- cbind(NFKB_smooth_all, treatment_id_all$treatment)

NFAT_with_ID_part3 <- subset(NFAT_with_ID_part2, NFAT_with_ID_part2$`treatment_id_all$treatment` %in% c("P1","P1+Ab"))
NFKB_with_ID_part3 <- subset(NFKB_with_ID_part2, NFKB_with_ID_part2$`treatment_id_all$treatment` %in% c("P1","P1+Ab"))


NFAT_with_ID <- rbind(NFAT_with_ID_part1, NFAT_with_ID_part3)
NFKB_with_ID <- rbind(NFKB_with_ID_part1, NFKB_with_ID_part3)

colnames(NFAT_with_ID) = colnames(NFKB_with_ID) = c(5*(1:48),"TREATMENT")

TREATMENT <- data.frame(table(NFAT_with_ID$TREATMENT))
TREATMENTs <- unique(TREATMENT$Var1)



#####
#################################### ___________________________________________________________________________________________
#### NFAT basal calculation ######## ________________________________________________________________________________________________________
#################################### ___________________________________________________________________________________________

NFAT_basal = data.frame()
Sub.min_data_NFAT = data.frame()
T4h.sd_max_NFAT = data.frame()

for (celli in 1:nrow(NFAT_with_ID)){ 
  
  dt.celli <- NFAT_with_ID[celli,1:48]
  NFAT_basal[celli,1] <- min(dt.celli) #apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean)
  
  Sub.min_data_NFAT[celli, 1:48]  <- NFAT_with_ID[celli,1:48] - NFAT_basal[celli,1]
  T4h.sd_max_NFAT[celli,1] <- apply(Sub.min_data_NFAT[celli,1:18],1, sd)
  T4h.sd_max_NFAT[celli,2] <- max(NFAT_with_ID[celli,1:18])/NFAT_basal[celli,1] ## within 4h
  T4h.sd_max_NFAT[celli,3] <- which.max(Sub.min_data_NFAT[celli,1:18])
  T4h.sd_max_NFAT[celli,4]  <- NFAT_with_ID$TREATMENT[celli] 
  T4h.sd_max_NFAT[celli,5] <- celli
  T4h.sd_max_NFAT[celli,6] <- apply((NFAT_with_ID[celli,1:48]/NFAT_basal[celli,1]-1),1, sum)
  T4h.sd_max_NFAT[celli,7] <- apply((NFAT_with_ID[celli,1:48]/NFAT_basal[celli,1]-1),1, max)
  
}

colnames(T4h.sd_max_NFAT) = c("transient.sd", "transient.max","transient.max.frame", "TREATMENT","cellID","AUC","transient.max.FC")
Sub.min_data_NFAT$TREATMENT <- NFAT_with_ID$TREATMENT

### NFKB basal calculation ########### ________________________________________________________________________________________________________

NFKB_basal = data.frame()
Sub.min_data_NFKB = data.frame()
T4h.sd_max_NFKB = data.frame()

for (celli in 1:nrow(NFKB_with_ID)){
  
  dt.celli <- NFKB_with_ID[celli,1:48]
  NFKB_basal[celli,1] <- apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean) #min(dt.celli)#
  
  Sub.min_data_NFKB[celli, 1:48] <- data.frame(NFKB_with_ID[celli,1:48] - NFKB_basal[celli,1])
  
  T4h.sd_max_NFKB[celli,1] <- apply(Sub.min_data_NFKB[celli,1:18],1, sd)
  T4h.sd_max_NFKB[celli,2] <- max(Sub.min_data_NFKB[celli,1:18])/NFKB_basal[celli,1]
  T4h.sd_max_NFKB[celli,3] <- which.max(Sub.min_data_NFKB[celli,1:18])
  T4h.sd_max_NFKB[celli,4] =  NFKB_with_ID$TREATMENT[celli]
  T4h.sd_max_NFKB[celli,5] <- celli
  T4h.sd_max_NFKB[celli,6] <- apply((NFKB_with_ID[celli,1:48]/NFKB_basal[celli,1]-1),1, sum)
  
}

colnames(T4h.sd_max_NFKB) = c("transient.sd", "transient.max","transient.max.frame","TREATMENT","cellID","AUC")

Sub.min_data_NFKB$TREATMENT <- NFKB_with_ID$TREATMENT

##### heatmaps (basal un) ######
library(pheatmap)
library(RColorBrewer)

##################
## Figure S5A,C ##
##################
# 
# pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/perturbation/allsc_P1Ab_P1_Ab_combined_ratio_2023.pdf",sep=""), width = 9,height = 3)
# # pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/perturbation/allsc_I0.6Ab_I0.6_Ab_combined_ratio_2023.pdf",sep=""), width = 9,height = 3)

  label_2 <- c("5","65","125","180","240")
    treatij <- c("P1+Ab", "Ab", "P1")
    # treatij <- c("I0.6+Ab", "Ab","I0.6")

  
  NFAT_treati = rbind(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]), subset(Sub.min_data_NFAT, TREATMENT == treatij[2]),subset(Sub.min_data_NFAT, TREATMENT == treatij[3]))
  NFKB_treati = rbind(subset(Sub.min_data_NFKB, TREATMENT == treatij[1]), subset(Sub.min_data_NFKB, TREATMENT == treatij[2]),subset(Sub.min_data_NFKB, TREATMENT == treatij[3]))
  
  # pht.A <-
  pheatmap(NFAT_treati[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),labels_col =c(0,120,240),
           main="NFAT under P1+Ab, Ab, P1",
           # main="NFAT under I0.6+Ab, Ab, I0.6",
           gaps_row=c(nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1])),
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2]))), #c("I0.3+Ab", "I0.3", "Ab")
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  # pht.B <-
  pheatmap(NFKB_treati[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),
           labels_col =c(0,120,240),
           main="NFKB under P1+Ab, Ab, P1",
           # main="NFKB under I0.6+Ab, Ab, I0.6",
           gaps_row=c(nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1])),
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2]))), #c("I0.3+Ab", "I0.3", "Ab")
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  dev.off()
  

##################
## Figure S5B,D ##
##################
### lineplot

library(reshape2)
library(dplyr)

treatijs <- list(c("I0.6+Ab", "Ab", "I0.6"), c("P1+Ab", "Ab", "P1"))
ipt <- c("Sub.min_data_NFAT", "Sub.min_data_NFKB")
title <- c("I0.6+Ab_Ab_I0.6","P1+Ab_Ab_P1")
TF <- c("NFAT", "NFKB")


ij=1
# ij=2
treatij <- treatijs[[ij]] # ij=1, Iono-nfat perturbation; # ij=2, PMA-nfkb perturbation

for (tf in 1:2){
  TF_treati <- subset(get(ipt[tf]), TREATMENT%in% treatij)
# #### only for ij=2 ##{
#   n = 0.40
#   NFAT_ranking0 <- subset(T4h.sd_max_NFAT, TREATMENT == "P1+Ab")
#   NFAT_ranking <- NFAT_ranking0[order(-NFAT_ranking0$AUC), ] 
#   rowid_NFAT_P1_Ab <- rownames(NFAT_ranking[1:(round(nrow(NFAT_ranking))*n),])
#   NFAT_ranking0 <- subset(T4h.sd_max_NFAT, TREATMENT == "Ab")
#   NFAT_ranking <- NFAT_ranking0[order(-NFAT_ranking0$AUC), ] 
#   rowid_NFAT_Ab <- rownames(NFAT_ranking[(round(nrow(NFAT_ranking))*(1-n)+1):nrow(NFAT_ranking),])
#   rowid_NFAT_P1 <- rownames(subset(T4h.sd_max_NFAT, TREATMENT == "P1"))
#   
#   TF_treati <- subset(TF_treati, rownames(TF_treati) %in% c(rowid_NFAT_P1_Ab, rowid_NFAT_Ab, rowid_NFAT_P1))
#   assign(paste("adj0.4_PAb_Ab_P_",TF[tf], sep=""),TF_treati)
# #### only for ij=2 ##}
  
  melt_data_TF <- reshape2::melt(TF_treati)
  melt_data_TF %>% group_by(TREATMENT, variable) %>% summarise(ave = mean(value), sd=sd(value)) -> TF_plt
  TF_plt$variable <- 5* as.numeric(TF_plt$variable)
  # TF_plt %>% group_by(TREATMENT) %>% summarise(max=max(ave), min=min(ave)) -> TF_plt_maxmin
  
  p<-
    ggplot(TF_plt, aes(x=variable, y=ave, group=TREATMENT, color = TREATMENT)) +
    geom_line(size=2)+ #rgb(49,140,231,maxColorValue=255)) +
    # geom_ribbon(aes(x=variable, ymin=ave-sd, ymax=ave+sd), fill=alpha(I(7/10)))+
    scale_x_continuous(breaks = seq(0, 240, 60))+
    labs(x = "Time(min)", y = "N/C ratio", title= paste(TF[tf],title[ij], sep=" under "))+
    themo_demo
  assign(paste("p",tf, sep="_"), p)
}

library(ggpubr)
# pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/perturbation/AVERAGE_lineplot_",title[ij],"_adjusted_0.4_2023.pdf",sep=""), width = 5,height = 5)
# pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/perturbation/AVERAGE_lineplot_",title[ij],"_2023.pdf",sep=""), width = 5,height = 5)
ggarrange(p_1,p_2, ncol =1,nrow =2)
dev.off()



###############
## Figure 2F ##
###############

################ analog modulation #########################

AUC.mean <- data.frame(treatment=0,NFAT.AUC=0,NFKB.AUC=0)

plot_all_transient.AB <- data.frame(cellID= paste("cell", rep(1:nrow(Sub.min_data_NFKB)), sep=""),
                                    TREATMENT = Sub.min_data_NFKB$TREATMENT,
                                    NFAT= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB= as.numeric(T4h.sd_max_NFKB$transient.max),
                                    NFAT.AUC= as.numeric(T4h.sd_max_NFAT$AUC),
                                    NFKB.AUC= as.numeric(T4h.sd_max_NFKB$AUC),
                                    NFAT.FC= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB.FC= as.numeric(T4h.sd_max_NFKB$transient.max))

##################  bootstarp #####################
library(bootstrap)

conf.int=function(x,sigma,conf.level=0.95) {
  mean<-mean(x)
  n<-length(x)
  alpha<-1-conf.level
  z=qnorm(1-alpha/2,mean=0,sd=1,lower.tail = T)
  c(mean-sigma*z/sqrt(n),mean+sigma*z/sqrt(n))
}

##
allR.mean <- data.frame()
allR.data <- data.frame()
allR.conf.pearson_cor <- data.frame()
cor.p.exp <-data.frame()
k=0

## for I0.6+Ab, I0.6
for (x in c(1,3,5)){ ## Ab, I0.6+Ab, P1+Ab
  if (x==5){## for adjusted 0.4 P1+Ab
     IDS <- "rowid_NFAT_P1_Ab"
     k=k+1
     treati = "P1+Ab"
     plot_x_AB <- subset(plot_all_transient.AB, rownames(plot_all_transient.AB) %in% get(IDS))
     cor.exp <- cor.test(plot_x_AB$NFAT.AUC, plot_x_AB$NFKB.AUC,method = "pearson")
     cor.p.exp[treati,1] <- cor.exp$p.value
  }else{    ## for I0.6+Ab, I0.6
    k=k+1
    treati = as.character(TREATMENTs[x])
    plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT == treati)
    cor.exp <- cor.test(plot_x_AB$NFAT.AUC, plot_x_AB$NFKB.AUC,method = "pearson")
    cor.p.exp[treati,1] <- cor.exp$p.value
  }
  
  # bootstrap cycle
  B <- 1000 #cycle
  n <- nrow(plot_x_AB) #size
  R <-  data.frame(treatment=0,cor.estimate=0,cor.p.value=0) #record correlation in each sampling process
  
  for (b in 1:B) {
    #randomly select the indices
    i <- sample(1:n, size = 0.5*n, replace = TRUE) # 1:n, replaceable
    plot_x_ABi <-  plot_x_AB[i,]
    aaa <- plot_x_ABi$NFAT.AUC
    bbb <- plot_x_ABi$NFKB.AUC
    cor.all <- cor.test(aaa, bbb,method = "pearson")
    R[b,1:3] <- c(treati,cor.all$estimate, cor.all$p.value)
  }
  sd.R=sd(as.numeric(R$cor.estimate))
  R.mean <- data.frame( treatment= treati, 
                        mean.R =mean(as.numeric(R$cor.estimate)),#mean pearson correlation
                        sd.R=sd(as.numeric(R$cor.estimate))) #sd pearson correlation
  
  allR.mean <- rbind(allR.mean, R.mean)
  allR.data <- rbind(allR.data, R)
  
  conf.pearson_cor <- data.frame(treatment= treati, L2=conf.int(as.numeric(R$cor.estimate),sd.R,0.95)[1], H2=conf.int(as.numeric(R$cor.estimate),sd.R,0.95)[2])  # conf.int(input_data, sample_sd, 0.95)
  allR.conf.pearson_cor <- rbind(allR.conf.pearson_cor, conf.pearson_cor)
}

allR.data$treatment = allR.mean$treatment
################

dt <- allR.mean
dt$L2 <- c(allR.conf.pearson_cor$L2)
dt$H2 <- c(allR.conf.pearson_cor$H2)
dt$p.value.exp <- cor.p.exp$V1



# pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/perturbation/Fig3F.Dyna_pearson_cor_adj0.4_2023.pdf",width = 6, height = 2)

library(ggplot2)

dt$treatment <- factor(dt$treatment, levels = c("I0.6+Ab","Ab","P1+Ab"), ordered = T)

  ggplot(data = dt, aes(x= treatment, y= mean.R, fill=p.value.exp)) +
  geom_bar(width = .5,stat="identity", position = 'dodge')+
  geom_errorbar(aes(ymin = L2, ymax = H2), width = 0.2, position = position_dodge(0.5))+
  scale_colour_gradient(low = "green", high = "red")+
  ylim(-0.04,0.3) +
  labs(title="pearson correlation",x="treatment",y="correlation")+
  # geom_signif(comparisons = compaired,step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo

dev.off()

