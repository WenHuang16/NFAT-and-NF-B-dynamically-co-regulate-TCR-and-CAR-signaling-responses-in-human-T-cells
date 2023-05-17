# ##__________________________________________________________________________________________________________________
# ## SINGLE-CELL DATA WERE USED AFTER DROPPING OUT OUTLIERS & SMOOTHING
# ##__________________________________________________________________________________________________________________
#
# library("R.matlab")
#
# ###### PMA & Iono ________________________________________________________________________________________________________
# ###### PI_matrix
# input_dir_PI_matrix <- "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/02.PMA and Ionomycin/metadata/"
#
# NFAT_EXCEL_smooth <- data.frame(readMat(paste(input_dir_PI_matrix, "NFAT_EXCEL2.mat", sep="")))
# NFKB_EXCEL_smooth <- data.frame(readMat(paste(input_dir_PI_matrix, "NFKB_EXCEL2.mat", sep="")))
#
# ID <-  read.csv(paste(input_dir_PI_matrix, "TREATMENT_ID.csv", sep=""), quote = "\"'",header = F)
# colnames(ID) = c("dir","file","date","pos","PMA","Iono","TREATMENT")
#
# PI_NFAT_with_ID <- cbind(NFAT_EXCEL_smooth[,1:96], ID[,7])
# PI_NFKB_with_ID <- cbind(NFKB_EXCEL_smooth[,1:96], ID[,7])
# colnames(PI_NFAT_with_ID) = colnames(PI_NFKB_with_ID) = c(5*(1:96),"TREATMENT")
#
#
# TREATMENT <- data.frame(table(ID$TREATMENT))
# colnames(TREATMENT) = c("ID","number")
# TREATMENT$ID = factor(TREATMENT$ID , levels = c(unique(ID$TREATMENT)), ordered = T)
#
# ###### I1.4P50 captured together with TCR antibody treatments
#
# NFAT_EXCEL_smooth_PI <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/I1.4P50_NFAT_smooth.mat", sep="")), TREATMENT="I1.4P50_add")
# NFKB_EXCEL_smooth_PI <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/I1.4P50_NFKB_smooth.mat", sep="")), TREATMENT="I1.4P50_add")
#
# NFAT_EXCEL_smooth_P <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/P50_NFAT_smooth.mat", sep="")), TREATMENT="I0P50_add")
# NFKB_EXCEL_smooth_P <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/P50_NFKB_smooth.mat", sep="")), TREATMENT="I0P50_add")
#
# NFAT_EXCEL_smooth_I <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/I1.4_NFAT_smooth.mat", sep="")), TREATMENT="I1.4P0_add")
# NFKB_EXCEL_smooth_I <- data.frame(readMat(paste(input_dir_PI_matrix, "PI_data_I1.4P50/I1.4_NFKB_smooth.mat", sep="")), TREATMENT="I1.4P0_add")
#
# PI_add_NFAT_with_ID  <- rbind(NFAT_EXCEL_smooth_PI[,c(1:96,109)], NFAT_EXCEL_smooth_P[,c(1:96,109)], NFAT_EXCEL_smooth_I[,c(1:96,109)])
# PI_add_NFKB_with_ID  <- rbind(NFKB_EXCEL_smooth_PI[,c(1:96,109)], NFKB_EXCEL_smooth_P[,c(1:96,109)], NFKB_EXCEL_smooth_I[,c(1:96,109)])
# colnames(PI_add_NFKB_with_ID) = colnames(PI_add_NFAT_with_ID) = colnames(PI_NFKB_with_ID)
#
#
#
# ###### CD3 & CD28 ________________________________________________________________________________________________________
# input_dir_TCR <- "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/01.TCR Antibody/metadata/"
#
# Ab_data_NFAT <- data.frame(readMat(paste(input_dir_TCR, "NFAT_EXCEL2.mat", sep="")))
# Ab_data_NFKB <-  data.frame(readMat(paste(input_dir_TCR, "NFKB_EXCEL2.mat", sep="")))
# colnames(Ab_data_NFAT)= colnames(Ab_data_NFKB) <- 5*(1:108)
#
# Ab_data_ID <-  read.csv(paste(input_dir_TCR, "Dragonfly_imaging_190911_22_ID_unique.csv", sep=""), quote="''",header=F)
# colnames(Ab_data_ID) <- c("PATH","CELL","TREATMENT")
#
# ID.frame <- data.frame(matrix(0,nrow(Ab_data_ID),5))
# colnames(ID.frame) <- c("PATH","CELL", "TREATMENT","DATE","POS")
#
# for (x in 1:nrow(Ab_data_ID))
# {
#   namei <- strsplit(as.character(Ab_data_ID$PATH[x]), "\\", fixed = T)
#   ID.frame[x,1:5] <- c(as.character(Ab_data_ID$PATH[x]),
#                        as.character(Ab_data_ID$CELL[x]),
#                        as.character(Ab_data_ID$TREATMENT[x]),
#                        namei[[1]][6],
#                        namei[[1]][8])
# }
#
# Ab_data_NFAT_ID <- cbind(Ab_data_NFAT[,1:97], ID.frame)
# Ab_data_NFKB_ID  <- cbind(Ab_data_NFKB[,1:97], ID.frame)
#
# #### Align imaging data from t1
# #### delete experiments containing t0
# Ab_data_NFAT_ID_with.T0 <- subset(Ab_data_NFAT_ID, Ab_data_NFAT_ID$DATE %in% c("2019-09-10", "2019-09-11", "2019-09-13", "2019-09-20"))
# Ab_data_NFKB_ID_with.T0 <- subset(Ab_data_NFKB_ID, Ab_data_NFKB_ID$DATE %in% c("2019-09-10", "2019-09-11", "2019-09-13", "2019-09-20"))
# NFAT_ID_del.T0 <- Ab_data_NFAT_ID_with.T0[,2:102]
# NFKB_ID_del.T0 <- Ab_data_NFKB_ID_with.T0[,2:102]
# colnames(NFAT_ID_del.T0) = colnames(NFKB_ID_del.T0) = c(5*(1:96),"PATH","CELL","TREATMENT","DATE","POS")
#
# ## combine all experiments without t0
# Ab_data_NFAT_ID_no.T0 <- subset(Ab_data_NFAT_ID, !Ab_data_NFAT_ID$DATE %in% c("2019-09-10", "2019-09-11", "2019-09-13", "2019-09-20"))
# Ab_data_NFKB_ID_no.T0 <- subset(Ab_data_NFKB_ID, !Ab_data_NFKB_ID$DATE %in% c("2019-09-10", "2019-09-11", "2019-09-13", "2019-09-20"))
#
# Ab_data_NFAT_ID_new.no.T0 <- rbind(Ab_data_NFAT_ID_no.T0[,c(1:96,98:102)],NFAT_ID_del.T0)
# Ab_data_NFKB_ID_new.no.T0 <- rbind(Ab_data_NFKB_ID_no.T0[,c(1:96,98:102)],NFKB_ID_del.T0)
#
#
# ###### combine TCR, PMA-Ionomycin matrix together ________________________________________________________________________________________________________
#
# All_data_NFKB_ID <- rbind(PI_NFKB_with_ID[,c(1:96,97)], PI_add_NFKB_with_ID[,c(1:96,97)], Ab_data_NFKB_ID_new.no.T0[,c(1:96,99)])
# All_data_NFAT_ID <- rbind(PI_NFAT_with_ID[,c(1:96,97)], PI_add_NFAT_with_ID[,c(1:96,97)], Ab_data_NFAT_ID_new.no.T0[,c(1:96,99)])
#
# ### save.image(paste(input_dir_TCR,"TCR_PI_T4h_combined_2023.Rdata", sep=""))

load("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_2023.Rdata")

library(ggplot2)
library(ggpubr)
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

# Ab_data_ID$date = Ab_data_ID$pos <- NA
#
# for (ri in 1:nrow(Ab_data_ID)){
#   strs <-  strsplit(as.character(Ab_data_ID$PATH[ri]),"[\\\\]|[^[:print:]]",fixed=FALSE)
#   Ab_data_ID$date[ri] <- strs[[1]][6]
#   Ab_data_ID$pos[ri] <- strs[[1]][7]
#   }
#
# table(Ab_data_ID[,c(4,3)])


### NFKB basal calculation ########### ________________________________________________________________________________________________________
# 5-90min: transient

NFKB_basal = data.frame()
Sub.min_data_NFKB = data.frame()
T4h.sd_max_NFKB = data.frame()

for (celli in 1:nrow(All_data_NFKB_ID)){

  dt.celli <- All_data_NFKB_ID[celli,1:96]
  NFKB_basal[celli,1] <- apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean)
  Sub.min_data_NFKB[celli, 1:48] <- data.frame(All_data_NFKB_ID[celli,1:48] - NFKB_basal[celli,1])

  T4h.sd_max_NFKB[celli,1] <- apply(Sub.min_data_NFKB[celli,1:18],1, sd)
  T4h.sd_max_NFKB[celli,2] <- max(Sub.min_data_NFKB[celli,1:18])/NFKB_basal[celli,1]
  T4h.sd_max_NFKB[celli,3] <- which.max(Sub.min_data_NFKB[celli,1:18])
  T4h.sd_max_NFKB[celli,4] =  All_data_NFKB_ID$TREATMENT[celli]
  T4h.sd_max_NFKB[celli,5] <- celli
  T4h.sd_max_NFKB[celli,6] <- apply((All_data_NFKB_ID[celli,1:48]/NFKB_basal[celli,1]-1),1, sum)

}

colnames(T4h.sd_max_NFKB) = c("transient.sd", "transient.max","transient.max.frame","TREATMENT","cellID","AUC")

########## add column to record single cell activation/not ######
plot_all_AB<- data.frame(cbind(cellID= paste("cell", rep(1:nrow(All_data_NFKB_ID)), sep=""), TREATMENT =as.character(All_data_NFKB_ID$TREATMENT)))
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

for (celli in 1:nrow(All_data_NFAT_ID)){

  dt.celli <- All_data_NFAT_ID[celli,1:96]
  NFAT_basal[celli,1] <- min(dt.celli) #apply(dt.celli[order(dt.celli,decreasing=F)[1:5]],1,mean)

  Sub.min_data_NFAT[celli, 1:96]  <- All_data_NFAT_ID[celli,1:96] - NFAT_basal[celli,1]
  T4h.sd_max_NFAT[celli,1] <- apply(Sub.min_data_NFAT[celli,1:18],1, sd)
  T4h.sd_max_NFAT[celli,2] <- max(All_data_NFAT_ID[celli,1:18])/NFAT_basal[celli,1] ## within 4h
  T4h.sd_max_NFAT[celli,3] <- which.max(Sub.min_data_NFAT[celli,1:18])
  T4h.sd_max_NFAT[celli,4]  <- All_data_NFAT_ID$TREATMENT[celli]
  T4h.sd_max_NFAT[celli,5] <- celli
  T4h.sd_max_NFAT[celli,6] <- apply((All_data_NFAT_ID[celli,1:48]/NFAT_basal[celli,1]-1),1, sum)
  T4h.sd_max_NFAT[celli,7] <- apply((All_data_NFAT_ID[celli,1:48]/NFAT_basal[celli,1]-1),1, max)

}

colnames(T4h.sd_max_NFAT) = c("transient.sd", "transient.max","transient.max.frame", "TREATMENT","cellID","AUC","transient.max.FC")
# plot(density(T4h.sd_max_NFAT$transient.max))


########## add column to record single cell activation/not ######
NFAT_max.thr.new  = max(subset(T4h.sd_max_NFAT, T4h.sd_max_NFAT$TREATMENT == "I0P0")[,6])
T4h.sd_max_NFAT$group <- ifelse((T4h.sd_max_NFAT$AUC > NFAT_max.thr.new),"response","non")

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


treat_order <- c("0.5CC", "CC", "2CC", "4CC",
                 "I0P0","I0.25P0","I0.5P0","I1P0","I1.4P0","I2.8P0",
                 "I0P5","I0.25P5","I0.5P5","I1P5","I1.4P5","I2.8P5",
                 "I0P25","I0.25P25","I0.5P25","I1P25","I1.4P25","I2.8P25",
                 "I0P50","I0.25P50","I0.5P50","I1P50","I1.4P50","I2.8P50",
                 "I0P75","I0.25P75","I0.5P75","I1P75","I1.4P75","I2.8P75",
                 "I0P125","I0.25P125","I0.5P125","I1P125","I1.4P125","I2.8P125",
                 "I0P250","I0.25P250","I0.5P250","I1P250","I1.4P250","I2.8P250",
                 "I1.4P50_add", "I0P50_add", "I1.4P0_add")


summary_abcd.transient <- data.frame()
for (treati in 1:length(treat_order)){

  plot_x_AB <- subset(plot_all_AB, plot_all_AB$TREATMENT == treat_order[treati])
  sum <- data.frame(table(plot_x_AB$transient.Color))
  sum$ratio <- sum$Freq/sum(sum$Freq)*100
  sum$TREATMENT <- treat_order[treati]
  summary_abcd.transient <- rbind(summary_abcd.transient, sum)

}


####################### TCR  activation  ##############################################

summary_abcd.TCR.all <- subset(summary_abcd.transient, summary_abcd.transient$TREATMENT %in% treat_order[1:4])
cols <- c("a"= "RoyalBlue", "b"="Green4","c"= "Tomato","d"= "Gray")
library(data.table)
library(reshape2)

TCR.plot.melt <- melt(summary_abcd.TCR.all[,c(1,3,4)])
TCR.plot.melt$TREATMENT <- factor(TCR.plot.melt$TREATMENT, levels = treat_order[1:4], ordered=T, )

###############
## Figure 1E ##
###############

ggplot(TCR.plot.melt) +
  geom_bar(aes(TREATMENT, value, fill = Var1), stat="identity", position="stack", alpha = 0.7, width=.5) +
  scale_fill_manual(values= cols) +
  ylim(0,100)+
  geom_hline(yintercept = 0)+
  labs(title = "")+
  themo_demo



####################### PI  activation  ##############################################

treat_order2 <- c("I0P0","I0.25P0","I0.5P0","I1P0",
                  "I0P5","I0.25P5","I0.5P5","I1P5",
                  "I0P25","I0.25P25","I0.5P25","I1P25",
                  "I0P50","I0.25P50","I0.5P50","I1P50",
                  "I1.4P50_add", "I0P50_add", "I1.4P0_add")

summary_abcd.PI.all <- subset(summary_abcd.transient, summary_abcd.transient$TREATMENT %in% treat_order2)

summary_abcd.PI.transient <- data.frame(TREATMENT = rep(treat_order2,4),
                                        Var1 = c(rep("a",length(treat_order2)),rep("b",length(treat_order2)),
                                                 rep("c",length(treat_order2) ),rep("d",length(treat_order2))))


for (nr in 1:nrow(summary_abcd.PI.transient)){
  nrxx <- subset(summary_abcd.PI.all,
                 summary_abcd.PI.all$TREATMENT == summary_abcd.PI.transient$TREATMENT[nr] &
                   summary_abcd.PI.all$Var1 == summary_abcd.PI.transient$Var1[nr])[,2:3]
  if (is.na(nrxx[1,1])){
    nrxx = data.frame(matrix(0,1,2))
    summary_abcd.PI.transient[nr,3:4] <- nrxx
  }else{
    summary_abcd.PI.transient[nr,3:4] <- nrxx
  }
}

colnames(summary_abcd.PI.transient) <- c("TREATMENT", "Var1", "Freq", "ratio")
summary_abcd.PI.transient$TREATMENT <- factor(summary_abcd.PI.transient$TREATMENT, levels =treat_order2,ordered = T)


######################## plot_all_transient.AB recorded activation strength/sd for dotplot ########################


plot_all_transient.AB <- data.frame(cellID= paste("cell", rep(1:nrow(All_data_NFKB_ID)), sep=""),
                                    TREATMENT = All_data_NFAT_ID$TREATMENT,
                                    NFAT= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB= as.numeric(T4h.sd_max_NFKB$transient.max),
                                    NFAT.AUC= as.numeric(T4h.sd_max_NFAT$AUC),
                                    NFKB.AUC= as.numeric(T4h.sd_max_NFKB$AUC),
                                    NFAT.FC= as.numeric(T4h.sd_max_NFAT$transient.max),
                                    NFKB.FC= as.numeric(T4h.sd_max_NFKB$transient.max),
                                    Color= plot_all_AB$transient.Color)

AUC.mean <- data.frame(NFAT.AUC=0,NFKB.AUC=0)
AUC.pos <- data.frame(NFAT.AUC=0,NFKB.AUC=0)
all.ratio <- data.frame(NFAT_pos.ratio=0, NFKB_pos.ratio=0)

treat_order3 <- c(treat_order[c(1:4)],treat_order2)

# table(PI_add_NFAT_with_ID$TREATMENT)

###########################
## Figure 1G, Figure S2F ##
###########################

for (treati in 1:length(treat_order3)){

  plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order3[treati])

  AUC.mean[treati,] <- data.frame(NFAT.AUC= mean(plot_x_AB$NFAT.AUC), NFKB.AUC= mean(plot_x_AB$NFKB.AUC))

  AUC.pos[treati,] <- data.frame(NFAT.AUC= mean(subset(plot_x_AB,Color %in% c("a", "b"))$NFAT.AUC),
                                 NFKB.AUC= mean(subset(plot_x_AB,Color %in% c("c", "b"))$NFKB.AUC))
  all.ratio[treati,] <- data.frame(NFAT_pos.ratio= nrow(subset(plot_x_AB,Color %in% c("a", "b")))/nrow(plot_x_AB),
                                   NFKB_pos.ratio=  nrow(subset(plot_x_AB, Color %in% c("c", "b")))/nrow(plot_x_AB))

  # cols <- c("a"= "RoyalBlue", "b"="SpringGreen3","c"= "Tomato","d"= "Gray")

  p <-
    ggplot(plot_x_AB) +
    geom_point(aes(NFKB.AUC, NFAT.AUC, color= Color, stroke = 0.8, alpha = 0.8))+
    scale_color_manual(values= cols) +
    xlim(0,100)+# xlim(-2,3)+#xlim(0,1)+
    ylim(0,148)+# ylim(0,4)+#ylim(0,6)+
    labs(title = treat_order[treati])+
    themo_demo

  assign(paste("p_",treat_order[treati],sep=""), p)

}
AUC.mean$treatment = AUC.pos$treatment = all.ratio$treatment = treat_order3

# write.csv(AUC.mean, "//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/For_Fig4.AUC_mean_TCR_PI_2023.csv", row.names = F)

###########################
  ## Figure S2B-D ##
###########################
  
  # pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/NFAT_NFKB_FigS2BCD_2023.pdf")
  
  ## response ratio
  all.ratio.TCR <- all.ratio[1:4,]
  all.ratio.TCR$treatment <- factor (all.ratio.TCR$treatment, level=c ("0.5CC","CC","2CC","4CC")) 
  
  ggplot(all.ratio.TCR, aes(x=treatment, y = NFKB_pos.ratio*100, fill = treatment)) +
   scale_fill_manual(name="Treatment", values=c("MistyRose1","PaleVioletRed2","VioletRed3","VioletRed4"))+
    geom_bar(width=0.3, stat="identity")+ ylim(0,100)+
    labs(x='treatment', y= 'percentage of cell with NFkB response')+
    themo_demo
  
  ggplot(all.ratio.TCR, aes(x=treatment, y = NFAT_pos.ratio*100, fill = treatment)) +
    scale_fill_manual(name="Treatment", values=c("LightSteelBlue","SteelBlue1","SteelBlue3","SteelBlue4"))+
    geom_bar(width=0.3, stat="identity")+ ylim(0,100)+
    labs(x='treatment', y= 'percentage of cell with NFAT response')+
    themo_demo
  
  ## peak height(FC)
  plot_all_transient.AB.TCR <- subset(plot_all_transient.AB,  TREATMENT %in% c("0.5CC","CC","2CC","4CC")) 
  
  ggplot(plot_all_transient.AB.TCR, aes(x=TREATMENT, y = NFKB.FC, fill = TREATMENT)) +
    geom_boxplot(aes(fill = TREATMENT), position = position_dodge(0.4),width=0.5, size = 0.25)+
    scale_fill_manual(name="Treatment", values=c("MistyRose1","PaleVioletRed2","VioletRed3","VioletRed4"))+
    labs(x='treatment', y= 'NFKB N/C ratio')+
    themo_demo
  
  ggplot(plot_all_transient.AB.TCR, aes(x=TREATMENT, y = NFAT.FC, fill = TREATMENT)) +
    geom_boxplot(aes(fill = TREATMENT), position = position_dodge(0.4),width=0.5, size = 0.25)+
    scale_fill_manual(name="Treatment", values=c("LightSteelBlue","SteelBlue1","SteelBlue3","SteelBlue4"))+
    labs(x='treatment', y= 'NFAT N/C ratio')+
    themo_demo
  
  ## % of NFKB activation cell with a 2nd peak
  
  transient.sd_max_NFKB = data.frame()
  stable.sd_max_NFKB = data.frame()
  
  for (celli in 1:nrow(All_data_NFKB_ID)){
    den_celli <- density(as.numeric(All_data_NFKB_ID[celli,1:48]))
    
    NFKB_basal[celli,1] <- min(All_data_NFKB_ID[celli,1:48])
    
    Sub.min_data_NFKB[celli, 1:96] <- data.frame(All_data_NFKB_ID[celli,1:48] - NFKB_basal[celli,1])
    
    transient.sd_max_NFKB[celli,1] <- apply(Sub.min_data_NFKB[celli,1:18],1, sd)
    stable.sd_max_NFKB[celli,1] <- apply(Sub.min_data_NFKB[celli,19:48],1, sd)
    
    transient.sd_max_NFKB[celli,2] <- max(Sub.min_data_NFKB[celli,1:48])/NFKB_basal[celli,1]
    stable.sd_max_NFKB[celli,2] <- max(Sub.min_data_NFKB[celli,19:48])/NFKB_basal[celli,1] 
    
    transient.sd_max_NFKB[celli,3] <- which.max(Sub.min_data_NFKB[celli,1:18])
    stable.sd_max_NFKB[celli,3] <- which.max(Sub.min_data_NFKB[celli,19:48]) +18
    
    transient.sd_max_NFKB[celli,4] = stable.sd_max_NFKB[celli,4] <- All_data_NFKB_ID$TREATMENT[celli]
    transient.sd_max_NFKB[celli,5] = stable.sd_max_NFKB[celli,5] <- celli
  }
  
  colnames(transient.sd_max_NFKB) = colnames(stable.sd_max_NFKB) = c("sd", "max", "max.frame", "TREATMENT","cellID")
  
  
  plot_all_AB.1st2nd <- data.frame(cbind(cellID= paste("cell", rep(1:nrow(All_data_NFKB_ID)), sep=""), TREATMENT = as.character(All_data_NFKB_ID$TREATMENT)))
  colnames(plot_all_AB.1st2nd) <- c("cellID","TREATMENT")
  
  NFKB_threshold = 0.1
  ######
  
  for (celli in 1:nrow(plot_all_AB.1st2nd)){
    if (transient.sd_max_NFKB$sd[celli] > 0.1 & transient.sd_max_NFKB$max.frame[celli] > 4){
      plot_all_AB.1st2nd$first.peak[celli] <- 1
    }else{
      plot_all_AB.1st2nd$first.peak[celli] <- 0
    }
    if ((stable.sd_max_NFKB$sd[celli] > 0.1 & stable.sd_max_NFKB$max.frame[celli] > 19 & stable.sd_max_NFKB$max.frame[celli] <58) &
        plot_all_AB.1st2nd$first.peak[celli] == 1){
      plot_all_AB.1st2nd$secd.peak[celli] <- 1
    }else{
      plot_all_AB.1st2nd$secd.peak[celli] <- 0
    }
  }
  
  sec.ratio <- subset(plot_all_AB.1st2nd[,2:4], TREATMENT %in% c("0.5CC","CC","2CC","4CC")) 
  library(dplyr)
  sec.ratio.TCR <- data.frame(
  sec.ratio %>% group_by(TREATMENT) %>% summarize(num= n(), first = sum(first.peak), sec = sum(secd.peak))
  )
  
  sec.ratio.TCR$TREATMENT <- factor (sec.ratio.TCR$TREATMENT, level=c ("0.5CC","CC","2CC","4CC")) 
  
  ggplot(sec.ratio.TCR, aes(x=TREATMENT, y = sec/num*100, fill = TREATMENT)) +
    scale_fill_manual(name="Treatment", values=c("MistyRose1","PaleVioletRed2","VioletRed3","VioletRed4"))+
    geom_bar(width=0.3, stat="identity")+ ylim(0,100)+
    labs(x='treatment', y= '% of cell with a 2nd peak')+
    themo_demo
  
  ggplot(sec.ratio.TCR, aes(x=TREATMENT, y = sec/first*100, fill = TREATMENT)) +
    scale_fill_manual(name="Treatment", values=c("MistyRose1","PaleVioletRed2","VioletRed3","VioletRed4"))+
    geom_bar(width=0.3, stat="identity")+ ylim(0,100)+
    labs(x='treatment', y= '% of NFKB activation cell with a 2nd peak')+
    themo_demo
  
  dev.off()
##


######################################################################################################################################################
save.image("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_AUC_2023.Rdata")
     load ("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_AUC_2023.Rdata")
######################################################################################################################################################
library(ggplot2)  


##### average line plot #####
treat_order4 <- c(treat_order3, "I0P0")
for ( ti in 1:length(treat_order4)){

  plt.NFAT <- subset(All_data_NFAT_ID, TREATMENT == treat_order4[ti])[,1:48]
  plt.NFKB <- subset(All_data_NFKB_ID, TREATMENT == treat_order4[ti])[,1:48]
  
  plt.NFAT_average = data.frame(time = c(5*(1:48)),mean= apply(plt.NFAT,2,mean), sd= apply(plt.NFAT,2,sd))
  plt.NFKB_average = data.frame(time = c(5*(1:48)),mean= apply(plt.NFKB,2,mean), sd= apply(plt.NFKB,2,sd))
  
  p.A <- ggplot(plt.NFAT_average, aes(x=time, y=mean)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), fill="LightSkyBlue", alpha = .5, linewidth= 3) +
    geom_line(size=1,color=rgb(30,144,255,maxColorValue=255)) +
    labs(x = "Time(min)", y = "NFAT N/C ratio", title = paste(treat_order4[ti], sep=""))+
    ylim(0,6)+
    scale_x_continuous(limits = c(0,240), breaks = c(0, 60,120,180,240)) + 
    themo_demo
  
  p.B <- ggplot(plt.NFKB_average, aes(x=time, y=mean)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), fill= "#FFE4B5", alpha = .5, linewidth = 3) +
    geom_line(size=1,color=rgb(255,127,0,maxColorValue=255)) +
    labs(x = "Time(min)", y = "NFKB N/C ratio", title = paste(treat_order4[ti], sep=""))+
    ylim(0,2)+
    scale_x_continuous(limits = c(0,240), breaks = c(0, 60,120,180,240)) + 
    themo_demo 
  assign(paste("p.A.", ti, sep=""), p.A)
  assign(paste("p.B.", ti, sep=""), p.B)
  
}


library(Rmisc)
# pdf(paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Average_lineplot_TCR_P_I_PI_2023.pdf", sep=""), width = 8,height = 4)
multiplot(p.A.1, p.A.2, p.A.3, p.A.4,
          p.B.1, p.B.2, p.B.3, p.B.4,
          layout=matrix(c(1:8), nrow=2, byrow=TRUE))

# # pdf(paste("x", sep=""), width = 10,height = 14)
# 
# multiplot(p.A.5, p.A.6, p.A.7, p.A.8, 
#           p.A.9, p.A.10,p.A.11, p.A.12, 
#           p.A.13, p.A.14, p.A.15, p.A.16,        
#           p.A.17, p.A.18, p.A.19, p.A.20,
#           layout=matrix(c(1:16), nrow=4, byrow=TRUE))
# 
# multiplot(p.B.5, p.B.6, p.B.7, p.B.8, 
#           p.B.9, p.B.10,p.B.11, p.B.12, 
#           p.B.13, p.B.14, p.B.15, p.B.16,
#           p.B.17, p.B.18, p.B.19, p.B.20,
#           layout=matrix(c(1:16), nrow=4, byrow=TRUE))

###############
## Figure 2C ##
###############

multiplot(p.A.24, p.A.21, p.A.22, p.A.23, 
          p.B.24, p.B.21, p.B.22, p.B.23,
          layout=matrix(c(1:8), nrow=2, byrow=TRUE))

# dev.off()




################## Responding Ratio, AUC  #####################
ABratio.PI <- data.frame(matrix(0,length(treat_order2),3))
for (treati in 1:length(treat_order2)){
  plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order2[treati])
  ABratio.PI[treati,1] <-treat_order2[treati]
  ABratio.PI[treati,2] <- nrow(subset(plot_x_AB, plot_x_AB$Color %in% c("a","b")))/nrow(plot_x_AB) #NFAT
  ABratio.PI[treati,3] <- nrow(subset(plot_x_AB, plot_x_AB$Color %in% c("c","b")))/nrow(plot_x_AB)
  
}

heatmap.A.ratio.PI = heatmap.B.ratio.PI <-  data.frame(matrix(0,4,4))
colnames(heatmap.A.ratio.PI) = colnames(heatmap.B.ratio.PI) <- c(0,0.25, 0.5, 1)
rownames(heatmap.A.ratio.PI)= rownames(heatmap.B.ratio.PI) <- c(0,5,25,50)
heatmap.A.ratio.PI[1,] <- ABratio.PI[1:4,2]
heatmap.A.ratio.PI[2,] <- ABratio.PI[5:8,2]
heatmap.A.ratio.PI[3,] <- ABratio.PI[9:12,2]
heatmap.A.ratio.PI[4,] <- ABratio.PI[13:16,2]
heatmap.B.ratio.PI[1,] <- ABratio.PI[1:4,3]
heatmap.B.ratio.PI[2,] <- ABratio.PI[5:8,3]
heatmap.B.ratio.PI[3,] <- ABratio.PI[9:12,3]
heatmap.B.ratio.PI[4,] <- ABratio.PI[13:16,3]


library(pheatmap)
library(RColorBrewer)

# pdf(paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/heatmap_PI.matrix_responding_ratio_2023.pdf", sep=""), width = 3,height = 3)

###############
## Figure 2D ##
###############

pheatmap(heatmap.A.ratio.PI, main="NFAT responding ratio",#annotation_row="PMA", annotation_col = "Ionomycin",
         cluster_cols=F, cluster_rows=F,height=0.2, width=0.2,
         color = colorRampPalette(brewer.pal(7,"YlGnBu"))(100))

pheatmap(heatmap.B.ratio.PI, main="NFKB responding ratio",#annotation_row="PMA", annotation_col = "Ionomycin",
         cluster_cols=F, cluster_rows=F,height=0.2, width=0.2,
         color = colorRampPalette(brewer.pal(7,"YlOrBr"))(100))

# dev.off()



################## Responding Ratio, AUC  #####################

for (treati in 1:length(treat_order3)){
  plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order3[treati])

  pA <- 
    ggplot(plot_x_AB, aes(x = NFAT.AUC)) + ggtitle(treat_order3[treati])+
    geom_density(alpha = 0.1, fill = "skyblue")+ geom_vline(xintercept=30,linetype = "dashed", colour = "grey" )+
    xlim(0,200)+themo_demo
  
  pB <- 
    ggplot(plot_x_AB, aes(x = NFKB.AUC)) + ggtitle(treat_order3[treati])+                        
    geom_density(alpha = 0.1, fill = "orange")+ geom_vline(xintercept=20,linetype = "dashed", colour = "grey")+
    xlim(0,75)+themo_demo
  
  assign(paste("pA",treat_order3[treati],sep="_"), pA)
  assign(paste("pB",treat_order3[treati],sep="_"), pB)

}



library(Rmisc)

##################
## Figure S4D-E ##
##################

# pdf(paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/density_PI.matrix_AUC_2023.pdf", sep=""), width = 8,height = 8)
multiplot(pA_I0P0,pA_I0.25P0,pA_I0.5P0,pA_I1P0,
          pA_I0P5,pA_I0.25P5,pA_I0.5P5,pA_I1P5,
          pA_I0P25,pA_I0.25P25,pA_I0.5P25,pA_I1P25,
          pA_I0P50,pA_I0.25P50,pA_I0.5P50,pA_I1P50,
          layout=matrix(c(1:16), nrow=4, byrow=TRUE))

multiplot(pB_I0P0,pB_I0.25P0,pB_I0.5P0,pB_I1P0,
          pB_I0P5,pB_I0.25P5,pB_I0.5P5,pB_I1P5,
          pB_I0P25,pB_I0.25P25,pB_I0.5P25,pB_I1P25,
          pB_I0P50,pB_I0.25P50,pB_I0.5P50,pB_I1P50,
          layout=matrix(c(1:16), nrow=4, byrow=TRUE))
dev.off()

##################  bootstarp #####################

library(bootstrap)

conf.int=function(x,sigma,conf.level=0.95) {
  mean<-mean(x)
  n<-length(x)
  alpha<-1-conf.level
  z=qnorm(1-alpha/2,mean=0,sd=1,lower.tail = T)
  c(mean-sigma*z/sqrt(n),mean+sigma*z/sqrt(n))
}

allR.mean <- data.frame()
allR.data <- data.frame()
allR.conf.pearson_cor <- data.frame()
cor.p.exp <-data.frame()
k=0
for (x in c(1:4,47:49)){ #
  k=k+1
  treati = treat_order[x]
  plot_x_AB <- subset(plot_all_transient.AB, plot_all_transient.AB$TREATMENT %in% treat_order[x])
  cor.exp <- cor.test(plot_x_AB$NFAT.AUC, plot_x_AB$NFKB.AUC,method = "pearson")
  cor.p.exp[treati,1] <- cor.exp$p.value

  ##bootstrap cycle
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
                        mean.KK=mean(as.numeric(R$V4)), #mean 
                        sd.KK =sd(as.numeric(R$V4)))  #sd 
  
  
  allR.mean <- rbind(allR.mean, R.mean)
  allR.data <- rbind(allR.data, R)
  
  conf.pearson_cor <- data.frame(treatment= treati, L2=conf.int(as.numeric(R$cor.estimate),sd.R,0.95)[1], H2=conf.int(as.numeric(R$cor.estimate),sd.R,0.95)[2])  # conf.int(input_data, sample_sd, 0.95)
  allR.conf.pearson_cor <- rbind(allR.conf.pearson_cor, conf.pearson_cor)

}

allR.data$treatment = allR.mean$treatment  #= allR.conf$treatment <- factor(allR.mean$treatment,levels = c("0.5CC","CC","2CC","4CC"))
################


dt <- allR.mean
dt$L2 <- c(allR.conf.pearson_cor$L2)
dt$H2 <- c(allR.conf.pearson_cor$H2)
dt$p.value.exp <- cor.p.exp$V1

# write.csv(dt,"//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Fig1H.Dyna_pearson_correlation_TCR_PI_data_2023.csv",row.names = F)


##################
## Figure 1F&1H ##
##################

# pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Fig1H.Dyna_pearson_correlation_TCR_doses_2023.pdf",width = 6, height = 2)
dt2 <- dt[c(1,2,4,3),]

library(ggplot2)
ggplot(data = dt2, aes(x= treatment, y= mean.R, fill=p.value.exp)) +
  geom_bar(width = .5,stat="identity", position = 'dodge')+
  geom_errorbar(aes(ymin = L2, ymax = H2), width = 0.2, position = position_dodge(0.5))+
  # ylim(0,0.6) +
  labs(title="pearson correlation",x="treatment",y="correlation")+
  # geom_signif(comparisons = compaired,step_increase = 0.1, map_signif_level = F,test = t.test)+
  themo_demo

ggplot(data = dt2, aes(x= treatment, y= mean.KK)) +
  geom_bar(width = .5,stat="identity", position = 'dodge')+
  geom_errorbar(aes(ymin =mean.KK-sd.KK, ymax =mean.KK+sd.KK), width = 0.2, position = position_dodge(0.5))+
  geom_hline(yintercept=1,linetype = "dashed", colour = "grey" )+
  ylim(0,1.5) +
  labs(title="active/passive",x="treatment",y="Measured / Random")+
  themo_demo

# dev.off()

# pdf("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Fig2E.Dyna_pearson_correlation_TCR_doses_2023.pdf",width = 6, height = 2)
dt3 <- dt[c(2,5),]

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

# dev.off()

