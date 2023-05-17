## CC
load("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_AUC_2023.Rdata")
Sub.min_data_NFAT$TREATMENT = Sub.min_data_NFKB$TREATMENT <- All_data_NFAT_ID$TREATMENT
Sub.min_data_NFAT.TCR <- subset(Sub.min_data_NFAT, TREATMENT =="CC")
Sub.min_data_NFKB.TCR <- subset(Sub.min_data_NFKB, TREATMENT =="CC")


##__________________________________________________________________________________________________________________
## SINGLE-CELL DATA WERE USED AFTER DROPPING OUT OUTLIERS & SMOOTHING
##__________________________________________________________________________________________________________________

library("R.matlab")
library(ggplot2)
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
###### 4 CARs ________________________________________________________________________________________________________
input_dir <- "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/04.Jurkat CAR-T/metadata/"

NFAT_EXCEL_smooth <- data.frame(readMat(paste(input_dir, "NFAT_EXCEL_smooth_2023.mat", sep="")))
NFKB_EXCEL_smooth <- data.frame(readMat(paste(input_dir, "NFKB_EXCEL_smooth_2023.mat", sep="")))

rawID <- data.frame(t(data.frame(readMat(paste(input_dir,"fls_all_data_2023.mat" ,sep="")))))
rawID$date <- NA
rawID$posID <- NA
rawID$CARID <- NA

for (ri in 1:nrow(rawID)){
  posnum <-  strsplit(as.character(rawID$folder[ri]),"pos")[[1]][2]
  rawID$posID[ri] <- strsplit(as.character(rawID$folder[ri]),"[\\\\]|[^[:print:]]",fixed=FALSE)[[1]][6]
  rawID$date[ri] <- strsplit(as.character(rawID$folder[ri]),"[\\\\]|[^[:print:]]",fixed=FALSE)[[1]][4]
  if (posnum %in% c(1:5)){
    rawID$CARID[ri] <- "CD19-recy-41BB"
  }else if(posnum %in% c(6:10)){
    rawID$CARID[ri] <- "CD19-WT-41BB"
  }else if(posnum %in% c(11:15)){
    rawID$CARID[ri] <- "CD19-CD28"
  } else if(posnum %in% c(16:20)){
    rawID$CARID[ri] <- "CD19-1st"
  }
}

# write.csv(table(rawID[,c(3,8)]),"//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/04.Jurkat CAR-T/metadata/table_CAR_date_pos.csv")

CAR_NFAT_with_ID <- cbind(NFAT_EXCEL_smooth[,1:96], rawID$CARID)
CAR_NFKB_with_ID <- cbind(NFKB_EXCEL_smooth[,1:96], rawID$CARID)
colnames(CAR_NFAT_with_ID) = colnames(CAR_NFKB_with_ID) = c(5*(1:96),"TREATMENT")

TREATMENT <- data.frame(table(rawID$CARID))
colnames(TREATMENT) = c("ID","number")
TREATMENT$ID = factor(TREATMENT$ID , levels = c("CD19-1st","CD19-CD28","CD19-WT-41BB","CD19-recy-41BB"), ordered = T)



### NFKB basal calculation ########### ________________________________________________________________________________________________________

NFKB_basal = data.frame()
Sub.min_data_NFKB = data.frame()
T4h.sd_max_NFKB = data.frame()

for (celli in 1:nrow(CAR_NFKB_with_ID )){
  dt.celli <- CAR_NFKB_with_ID[celli,1:96]
  NFKB_basal[celli,1] <- apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean)
  
  Sub.min_data_NFKB[celli, 1:48] <- data.frame(CAR_NFKB_with_ID [celli,1:48] - NFKB_basal[celli,1])
  
  T4h.sd_max_NFKB[celli,1] <- apply(Sub.min_data_NFKB[celli,1:18],1, sd)
  T4h.sd_max_NFKB[celli,2] <- max(Sub.min_data_NFKB[celli,1:18])/NFKB_basal[celli,1]
  T4h.sd_max_NFKB[celli,3] <- which.max(Sub.min_data_NFKB[celli,1:18])
  T4h.sd_max_NFKB[celli,4] =  CAR_NFKB_with_ID $TREATMENT[celli]
  T4h.sd_max_NFKB[celli,5] <- celli
  T4h.sd_max_NFKB[celli,6] <-  apply((CAR_NFKB_with_ID[celli,1:48]/NFKB_basal[celli,1]-1),1, sum)
  
}

colnames(T4h.sd_max_NFKB) = c("transient.sd", "transient.max","transient.max.frame","TREATMENT","cellID","AUC")


########## add response column to record single cell activation/not ######
plot_all_AB<- data.frame(cellID= paste("cell", rep(1:nrow(CAR_NFKB_with_ID )), sep=""), TREATMENT =CAR_NFKB_with_ID$TREATMENT)
colnames(plot_all_AB) <- c("cellID","TREATMENT")

NFKB_threshold = 0.1
######

for (celli in 1:nrow(plot_all_AB)){
  if (T4h.sd_max_NFKB$transient.sd[celli] > 0.1 & T4h.sd_max_NFKB$transient.max.frame[celli] > 4){
    plot_all_AB$first.peak[celli] <- 1
  }else{
    plot_all_AB$first.peak[celli] <- 0
  }
}



#################################### ___________________________________________________________________________________________
#### NFAT basal calculation ######## ________________________________________________________________________________________________________
#################################### ___________________________________________________________________________________________

NFAT_basal = data.frame()
Sub.min_data_NFAT = data.frame()
T4h.sd_max_NFAT = data.frame()

for (celli in 1:nrow(CAR_NFAT_with_ID)){ 
  
  dt.celli <- CAR_NFAT_with_ID[celli,1:96]
  NFAT_basal[celli,1] <- apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean)


  Sub.min_data_NFAT[celli, 1:48]  <- CAR_NFAT_with_ID[celli,1:48] - NFAT_basal[celli,1]
  T4h.sd_max_NFAT[celli,1] <- apply(Sub.min_data_NFAT[celli,1:18],1, sd)
  T4h.sd_max_NFAT[celli,2] <- max(CAR_NFAT_with_ID[celli,1:18])/NFAT_basal[celli,1] ## within 4h
  T4h.sd_max_NFAT[celli,3] <- which.max(Sub.min_data_NFAT[celli,1:18])
  T4h.sd_max_NFAT[celli,4]  <- CAR_NFAT_with_ID$TREATMENT[celli]
  T4h.sd_max_NFAT[celli,5] <- celli
  T4h.sd_max_NFAT[celli,6] <- apply((CAR_NFAT_with_ID[celli,1:48]/NFAT_basal[celli,1]-1),1, sum)
}

colnames(T4h.sd_max_NFAT) = c("transient.sd", "transient.max","transient.max.frame", "TREATMENT","cellID","AUC")

###### transient activation ratio ##############################
NFAT_max.thr.new  = 2.004285 # max(subset(T4h.sd_max_NFAT, T4h.sd_max_NFAT$TREATMENT == "I0P0")[,2]) # medium basal

x <- subset(T4h.sd_max_NFAT, T4h.sd_max_NFAT$transient.max > NFAT_max.thr.new)
T4h.sd_max_NFAT$group <- ifelse((T4h.sd_max_NFAT$transient.max > NFAT_max.thr.new & T4h.sd_max_NFAT$transient.max.frame > 1) ,"response","non")


T4h.positive_NFAT <-  merge(data.frame(table(subset(T4h.sd_max_NFAT,group=="response")$TREATMENT)), data.frame(table(CAR_NFAT_with_ID$TREATMENT)), by="Var1",all = T )
T4h.positive_NFAT$percent <- T4h.positive_NFAT$Freq.x/T4h.positive_NFAT$Freq.y

# ###################### peak location
for (celli in 1:nrow(plot_all_AB)){
  if (T4h.sd_max_NFAT$group[celli] == "response"){
    plot_all_AB$NFAT.transient[celli] <-1
  }else{
    plot_all_AB$NFAT.transient[celli] <-0
  }
}


######################## plot_all_AB recorded activation(1) or not(0) ##########################

for (celli in 1:nrow(plot_all_AB)){
  
  x = plot_all_AB$first.peak[celli]
  y = plot_all_AB$NFAT.transient[celli]
  
  if (y ==1 & x ==0){
    plot_all_AB$transient.Color[celli] = "a"}
  else if(y ==1 & x ==1){
    plot_all_AB$transient.Color[celli] = "b"}
  else if (y ==0 & x ==1){
    plot_all_AB$transient.Color[celli] = "c"}
  else {
    plot_all_AB$transient.Color[celli] = "d"}
} 

plot_all_AB$NFAT.AUC <- T4h.sd_max_NFAT$AUC
plot_all_AB$NFKB.AUC <- T4h.sd_max_NFKB$AUC

treat_order <- c("CD19-1st", "CD19-CD28", "CD19-WT-41BB", "CD19-recy-41BB")

summary_abcd.transient <- data.frame()
for (treati in 1:length(treat_order)){
  
  plot_x_AB <- subset(plot_all_AB, plot_all_AB$TREATMENT == treat_order[treati])
  sum <- data.frame(table(plot_x_AB$transient.Color))
  sum$ratio <- sum$Freq/sum(sum$Freq)*100
  sum$TREATMENT <- treat_order[treati]
  summary_abcd.transient <- rbind(summary_abcd.transient, sum)
  
}


####################### CAR ativation ############################################## 

summary_abcd.TCR.all <- subset(summary_abcd.transient, summary_abcd.transient$TREATMENT %in% treat_order[1:4])
summary_abcd.TCR.all <- rbind(summary_abcd.TCR.all, c("c",0,0,"CD19-1st"))

summary_abcd.TCR.all$Freq <- as.numeric(summary_abcd.TCR.all$Freq)
####
cols <- c("a"= "RoyalBlue", "b"="Green4","c"= "Tomato","d"= "Gray")
library(data.table)

TCR.plot.melt <- reshape2::melt(summary_abcd.TCR.all[,c(1,3,4)]) 
TCR.plot.melt$TREATMENT <- factor(TCR.plot.melt$TREATMENT, levels = treat_order[1:4], ordered=T, )


##################
## Figure 3F ##
##################
# pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3F_CAR_transient.a.b.c.d_percentage.pdf")

TCR.plot.melt$Var1 = factor(TCR.plot.melt$Var1,levels=c("a","b","c","d"))
p0<-
  ggplot(TCR.plot.melt) +
  geom_bar(aes(TREATMENT, as.numeric(ratio), fill =Var1 ), stat="identity", position="stack", alpha = 0.7, width=.5) +
  scale_fill_manual(values= cols) +
  ylim(0,101)+
  geom_hline(yintercept = 0)+
  labs(title = "")+
  themo_demo

# dev.off()


plot_all_transient.AB <- data.frame(cellID= paste("cell", rep(1:nrow(CAR_NFAT_with_ID)), sep=""),
                                    TREATMENT = CAR_NFAT_with_ID$TREATMENT,
                                    NFAT= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB= as.numeric(T4h.sd_max_NFKB$transient.max),
                                    NFAT.AUC= as.numeric(T4h.sd_max_NFAT$AUC),
                                    NFKB.AUC= as.numeric(T4h.sd_max_NFKB$AUC),
                                    NFAT.FC= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB.FC= as.numeric(T4h.sd_max_NFKB$transient.max),
                                    Color= plot_all_AB$transient.Color)


AUC.mean.CAR <- data.frame(NFAT.AUC=0,NFKB.AUC=0)

for (treati in 1:4){#length(treat_order)){

  plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order[treati])
  AUC.mean.CAR[treati,] <- data.frame( NFAT.AUC= mean(plot_x_AB$NFAT.AUC), NFKB.AUC= mean(plot_x_AB$NFKB.AUC))

  cols <- c("a"= "RoyalBlue", "b"="SpringGreen3","c"= "Tomato","d"= "Gray")
  
  p <-
  ggplot(plot_x_AB) +
    geom_point(aes(NFKB.AUC, NFAT.AUC, color= Color, stroke = 0.8, alpha = 0.8))+
    scale_color_manual(values= cols) +
    xlim(0,100)+
    ylim(0,500)+
    labs(title = treat_order[treati])+
    themo_demo
  
    assign(paste("p_",treati,sep=""), p)

}


AUC.mean.CAR$treatment <- treat_order
# write.csv(AUC.mean.CAR, "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/For.Fig4.AUC_mean_CAR_2023.csv", row.names = F)
  
  
################
## Figure S7A ##
################
library(ggpubr)
# pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/FigS7A_CAR_AUC_dotplot.pdf", width =10,height = 12)
ggarrange( p_1,p_2,p_3,p_4,ncol =2,nrow =2)
# dev.off()

#######################
## Figure 3D heatmap ##
#######################

# pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3D_3CARs_heatmap.pdf",sep=""), width = 9,height = 3)
library(pheatmap)  
library(RColorBrewer)  
library(pheatmap)  

  label_2 <- c("5","65","125","180","240")
  treatij <- c("CD19-1st", "CD19-CD28", "CD19-WT-41BB","CC")
 
  Sub.min_data_NFAT$TREATMENT = Sub.min_data_NFKB$TREATMENT <- CAR_NFAT_with_ID$TREATMENT
  NFAT_treati = rbind(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]), subset(Sub.min_data_NFAT, TREATMENT == treatij[2]),subset(Sub.min_data_NFAT, TREATMENT == treatij[3]),subset(Sub.min_data_NFAT, TREATMENT == treatij[4]))
  NFKB_treati = rbind(subset(Sub.min_data_NFKB, TREATMENT == treatij[1]), subset(Sub.min_data_NFKB, TREATMENT == treatij[2]),subset(Sub.min_data_NFKB, TREATMENT == treatij[3]),subset(Sub.min_data_NFKB, TREATMENT == treatij[4]))
  
  # pht.A <-
  pheatmap(NFAT_treati[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),labels_col =c(0,120,240),
           main="NFAT under CD19-1st, CD19-CD28, CD19-WT-41BB",
           gaps_row=c(nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1])),
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2])),
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[3]))), #c("I0.3+Ab", "I0.3", "Ab")
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  # pht.B <-
  pheatmap(NFKB_treati[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),
           labels_col =c(0,120,240),
           main="NFKB under CD19-1st, CD19-CD28, CD19-WT-41BB",
           gaps_row=c(nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1])),
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2])), #c("I0.3+Ab", "I0.3", "Ab")
                      nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[1]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[2]))+nrow(subset(Sub.min_data_NFAT, TREATMENT == treatij[3]))), #c("I0.3+Ab", "I0.3", "Ab")
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  dev.off()
  
  
#######################
## Figure 3E ##
#######################  

#### average lineplot
library(reshape2)
library(dplyr)

  
ipt.CAR <- c("Sub.min_data_NFAT", "Sub.min_data_NFKB")
ipt.TCR <- c("Sub.min_data_NFAT.TCR", "Sub.min_data_NFKB.TCR")
title <- "CD19-1st, CD19-CD28, CD19-WT-41BB, CC"
TF <- c("NFAT", "NFKB")
  
treatij <- c("CD19-1st", "CD19-CD28", "CD19-WT-41BB","CC")
  
for (tf in 1:2){
  CAR_i <- get(ipt.CAR[tf])
  CC_i <- get(ipt.TCR[tf])[,c(1:48,97)]
  colnames(CAR_i) = colnames(CC_i) = c(5*seq(1:48),"TREATMENT")
  TF_treati <- subset(rbind(CAR_i,CC_i),TREATMENT%in% treatij)

    melt_data_TF <- reshape2::melt(TF_treati)
    melt_data_TF %>% group_by(TREATMENT, variable) %>% summarise(ave = mean(value), sd=sd(value)) -> TF_plt
    TF_plt$variable <- 5*as.numeric(TF_plt$variable)
    TF_plt_maxmin <- TF_plt %>% group_by(TREATMENT) %>% summarise(max=max(ave), min=min(ave)) 
    
    p<- 
      ggplot(TF_plt, aes(x=variable, y=as.numeric(ave), group=TREATMENT, color = TREATMENT)) +
      geom_line(size=2)+ 
      scale_x_continuous(breaks = seq(0, 240, 60))+
      labs(x = "Time(min)", y = "N/C ratio", title= paste(TF[tf],title, sep=" under "))+
      themo_demo
    assign(paste("p",tf, sep="_"), p)
}

treat_CAR <- treatij[1:3]
for ( ti in 1:length(treat_CAR)){
  
  plt.NFAT <- subset(Sub.min_data_NFAT, TREATMENT == treat_CAR[ti])[,1:48]
  plt.NFKB <- subset(Sub.min_data_NFKB, TREATMENT == treat_CAR[ti])[,1:48]
  
  plt.NFAT_average = data.frame(time = c(5*(1:48)),mean= apply(plt.NFAT,2,mean), sd= apply(plt.NFAT,2,sd))
  plt.NFKB_average = data.frame(time = c(5*(1:48)),mean= apply(plt.NFKB,2,mean), sd= apply(plt.NFKB,2,sd))
  
  p.A <- ggplot(plt.NFAT_average, aes(x=time, y=mean)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), fill="LightSkyBlue", alpha = .5, linewidth= 3) +
    geom_line(size=1,color=rgb(30,144,255,maxColorValue=255)) +
    labs(x = "Time(min)", y = "NFAT N/C ratio", title = paste(treat_CAR[ti], sep=""))+
    ylim(-1,6)+
    scale_x_continuous(limits = c(0,240), breaks = c(0, 60,120,180,240)) + 
    themo_demo
  
  p.B <- ggplot(plt.NFKB_average, aes(x=time, y=mean)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), fill= "#FFE4B5", alpha = .5, linewidth = 3) +
    geom_line(size=1,color=rgb(255,127,0,maxColorValue=255)) +
    labs(x = "Time(min)", y = "NFKB N/C ratio", title = paste(treat_CAR[ti], sep=""))+
    ylim(-0.1,1)+
    scale_x_continuous(limits = c(0,240), breaks = c(0, 60,120,180,240)) + 
    themo_demo 
  assign(paste("p.A.", ti, sep=""), p.A)
  assign(paste("p.B.", ti, sep=""), p.B)
  
}

library(ggpubr)
# pdf(file = paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3E_AVERAGE_3CARS_CC.pdf",sep=""), width = 5,height = 5)
ggarrange(p.A.1,p.A.2,p.A.3,p.B.1,p.B.2,p.B.3, ncol =3,nrow =2)
ggarrange(p_1,p_2, ncol =1,nrow =2)
# dev.off()
  
  
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
  for (x in 1:4){ #
    k=k+1
    treati = treat_order[x]
    plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order[x])
    cor.exp <- cor.test(plot_x_AB$NFAT.AUC, plot_x_AB$NFKB.AUC,method = "pearson")
    cor.p.exp[treati,1] <- cor.exp$p.value
    
    #bootstrap cycle
    B <- 1000 #cycle
    n <- nrow(plot_x_AB) #size
    R <-  data.frame(treatment=0,cor.estimate=0,cor.p.value=0) #record correlation in each sampling process
    
    for (b in 1:B) {
      #randomly select the indices
      i <- sample(1:n, size = 0.5*n, replace = TRUE) # 1:n,replaceable
      plot_x_ABi <-  plot_x_AB[i,]
      aaa <- plot_x_ABi$NFAT.AUC
      bbb <- plot_x_ABi$NFKB.AUC
      
      cor.all <- cor.test(aaa, bbb,method = "pearson")
      R[b,1:3] <- c(treati,cor.all$estimate, cor.all$p.value)
      
      ## active VS. passive
      NFAT.P <- nrow(subset(plot_x_ABi, Color %in% c("a","b")))/(0.5*n)
      NFKB.P <- nrow(subset(plot_x_ABi, Color %in% c("c","b")))/(0.5*n)
      passive.DP <- NFAT.P*NFKB.P
      active.DP <- nrow(subset(plot_x_ABi, Color %in% c("b")))/(0.5*n)
      KK <- active.DP/passive.DP
      R[b,4] <- KK
      
    }

    sd.R=sd(as.numeric(R$cor.estimate))
    R.mean <- data.frame( treatment= treati, 
                          mean.R =mean(as.numeric(R$cor.estimate)),#mean pearson correlation
                          sd.R=sd(as.numeric(R$cor.estimate)), #sd pearson correlation
                          mean.KK=mean(na.omit(as.numeric(R$V4))), #mean 
                          sd.KK =sd(na.omit(as.numeric(R$V4))))  #sd 
    
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
  dt2 <- dt[c(1,2,3),]  
  
  # write.csv(dt,"//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3H.Dyna_pearson_correlation_CAR_data.csv",row.names = F)


  
##########  Figure 3 TCR & CAR #######################################################################################################################
pearson_correlation_TCR_PI <- read.csv("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Fig1H.Dyna_pearson_correlation_TCR_PI_data_2023.csv")
pearson_correlation_CAR <- read.csv("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3H.Dyna_pearson_correlation_CAR_data.csv")
  
dt3 <- rbind(pearson_correlation_TCR_PI[2,], pearson_correlation_CAR[1:3,])
  
pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig3/Fig3H.Dyna_pearson_correlation_TCR_CARs_2023.pdf",width = 6, height = 2)
  
  library(ggplot2)
  ggplot(data = dt3, aes(x= treatment, y= mean.R, fill=p.value.exp)) +
    geom_bar(width = .5,stat="identity", position = 'dodge')+
    geom_errorbar(aes(ymin = L2, ymax = H2), width = 0.2, position = position_dodge(0.5))+
    # ylim(0,0.6) + 
    labs(title="pearson correlation",x="treatment",y="correlation")+
    # geom_signif(comparisons = compaired,step_increase = 0.1, map_signif_level = F,test = t.test)+
    themo_demo
  
  ggplot(data = dt3, aes(x= treatment, y= mean.KK)) +
    geom_bar(width = .5,stat="identity", position = 'dodge')+
    geom_errorbar(aes(ymin =mean.KK-sd.KK, ymax =mean.KK+sd.KK), width = 0.2, position = position_dodge(0.5))+
    geom_hline(yintercept=1,linetype = "dashed", colour = "grey" )+
    ylim(0,1.5) +
    labs(title="active/passive",x="treatment",y="Measured / Random")+
    themo_demo
  
  dev.off()
  