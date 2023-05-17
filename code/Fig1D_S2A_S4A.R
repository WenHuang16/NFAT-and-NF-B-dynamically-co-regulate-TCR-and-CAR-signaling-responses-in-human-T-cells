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

load("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Imaging/01.TCR Antibody/metadata/TCR_PI_T4h_combined_2023.Rdata")

library("pheatmap")
library("RColorBrewer")

themo_demo= theme_bw()+
  theme(
    text = element_text(size=8),
    panel.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text=element_text(color='black'),
    plot.title = element_text(hjust = 0.5),
    title=element_text(size = 8),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8))


# # cell-treatment ID under anti-CD3&CD28 stimulation: "ID.frame"
# # final processed N/C ration: "Ab_data_NFKB_ID_new.no.T0", "Ab_data_NFAT_ID_new.no.T0"
table.Ab <- data.frame(table(ID.frame$DATE,ID.frame$TREATMENT ))
# 
# # cell-treatment ID under PMA, Iono, or P&I stimulation: "ID"
# # final processed N/C ration: "PI_add_NFKB_with_ID", "PI_add_NFAT_with_ID"
table.PI <- data.frame(table(ID$date, ID$TREATMENT))


########


########################################### ggplot2.two_y_axis ######################################################
##gtable + grid
library(ggplot2)
library(gtable)
library(grid) 

## define
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}




########################################### pLOT WITH SINGLE CELL ######################################################

library("ggplot2") 
library("Rmisc")
library("plyr")

tsAb<- unique(as.character(table.Ab$Var2))
tsPI <- unique(as.character(table.PI$Var2))
treatment <- c(tsAb, tsPI)

# Make plots. 


for (x in 5:length(treatment)){
  treati = treatment[x]  
  
  # pdf(file = paste("D:/LabNAS_Backup/HuangWen/00. 20221101_HW_paper_reviewed/allsc_",treati,"_NC_ratio.pdf",sep=""), width = 9,height = 3)
  label_2 <- c("5","65","125","180","240")
  

  
  NFAT_treati = subset(All_data_NFAT_ID, All_data_NFAT_ID$TREATMENT == treati)
  NFKB_treati = subset(All_data_NFKB_ID, All_data_NFKB_ID$TREATMENT == treati)
  
  for (i in 1:nrow(NFKB_treati)) {
    
    raw_d = data.frame(time= 5*(1:47),
                       NC_ratio_NFAT= t(NFAT_treati[i,1:47]),
                       NC_ratio_NFKB= t(NFKB_treati[i,1:47]))


    colnames(raw_d) <- c("time","NFAT_NC", "NFKB_NC")
    
    
    p1<- ggplot(raw_d, aes(x=time, y=raw_d$NFAT_NC)) +
      geom_line(size=2,color=rgb(30,144,255,maxColorValue=255)) +
      labs(x = "Time(min)", y = "N/C ratio",
           title = paste("NFAT ratio under ",treati,"_cell_",i,sep=""))+
      ylim(0,4)+scale_x_continuous(breaks = c(5,65,125,185,240), labels = label_2)+
      theme_bw()+
      theme(
        text = element_text(size=12),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color='black'),
        plot.title = element_text(hjust = 0.5),
        title=element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
    
    
    
    p2<-  ggplot(raw_d, aes(x=time, y=raw_d$NFKB_NC)) +
      geom_line(size=2,color=rgb(243,156,0,maxColorValue=255)) +
      labs(x = "Time(min)", y = "N/C ratio",
           title = paste("NFKB N/C ratio under ",treati,"_cell_",i,sep=""))+
      ylim(0,1.5)+scale_x_continuous(breaks = c(0,60,120,180,240), labels = label_2)+
      theme_bw()+
      theme(
        text = element_text(size=12),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color='black'),
        plot.title = element_text(hjust = 0.5),
        title=element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
    ###_____________________________________________________________________________________________
    
    p3 <- ggplot2.two_y_axis(p1, p2)
    
    multiplot(p1,p2,p3,layout=matrix(c(1,2), nrow=1, byrow=TRUE))
  }
  
  # dev.off() 
}



####################################### heatmap ###################################################################################
#### TCR

################
## Figure S2A ##
################
tsAb<- unique(as.character(table.Ab$Var2))
for (i in 1:4){

  pdf(paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/heatmap_TCR_", tsAb[i], "_2023.pdf", sep=""))
  
  pht.data_NFAT <- subset(All_data_NFAT_ID, All_data_NFAT_ID$TREATMENT == tsAb[i])
  # pht.A <-
    pheatmap(pht.data_NFAT[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),
           labels_col =c(0,120,240),
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  pht.data_NFKB <- subset(All_data_NFKB_ID, All_data_NFKB_ID$TREATMENT ==  tsAb[i])
  # pht.B <-
    pheatmap(pht.data_NFKB[,1:48], scale = "row", #clustering_method = "ward.D2",
             color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
             legend_breaks=seq(-4,4,2),
             labels_col =c(0,120,240),
             cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
             show_rownames = F, show_colnames = T)
  
  dev.off()
}

#### PI 
##################
## Figure S4A-B ##
##################

tsPI <- c("I0P0","I0.25P0","I0.5P0","I1P0",
          "I0P5","I0.25P5","I0.5P5","I1P5",
          "I0P25","I0.25P25","I0.5P25","I1P25",
          "I0P50","I0.25P50","I0.5P50","I1P50")
for (j in 1:length(tsPI)){

  pdf(paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/heatmap_PI.matrix_", tsPI[j], "_2023.pdf", sep=""))
  
    pht.data_NFAT <- subset(All_data_NFAT_ID, All_data_NFAT_ID$TREATMENT == tsPI[j])
  # pht.A <-
  pheatmap(pht.data_NFAT[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),
           labels_col =c(0,120,240),
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  pht.data_NFKB <- subset(All_data_NFKB_ID, All_data_NFKB_ID$TREATMENT ==  tsPI[j])
  # pht.B <-
  pheatmap(pht.data_NFKB[,1:48], scale = "row", #clustering_method = "ward.D2",
           color =  colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(100),
           legend_breaks=seq(-4,4,2),
           labels_col =c(0,120,240),
           cluster_cols = F, cluster_rows = F, cellwidth =1, cellheight = 0.1,cutree_rows = 2,
           show_rownames = F, show_colnames = T)
  
  dev.off()
}





