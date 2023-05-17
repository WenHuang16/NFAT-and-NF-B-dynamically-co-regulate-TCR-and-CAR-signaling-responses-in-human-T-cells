
library(pheatmap)
library(RColorBrewer)
# bk <- c(seq(0,1,by=0.01))
# color = c(colorRampPalette(colors = c("Snow","Deeppink4"))(length(bk)))
# 

setwd("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/FACS/JNAB-AB_or_BA/analysis_results/")

####
monoclone <- c("AB_1","AB_3","BA_4")
# i = 2
i = 3

pNFKB_ratio_8h <- read.csv(paste(monoclone[i],"_ratio_pNFKB_8h.csv",sep=""), header =F, sep = ",")
pNFAT_ratio_8h <- read.csv(paste(monoclone[i],"_ratio_pNFAT_8h.csv",sep=""), header =F, sep = ",")
DP_ratio_8h <-  read.csv(paste(monoclone[i],"_ratio_DP_8h.csv",sep=""), header =F, sep = ",")
mCherry_8h <- read.csv(paste(monoclone[i],"_ratio_mCherry_8h.csv",sep=""), header =F, sep = ",")

rownames(pNFKB_ratio_8h) = rownames(pNFAT_ratio_8h) = rownames(DP_ratio_8h) = rownames(mCherry_8h) = c("2.8","1.4", "1","0.5","0.25", "0")
colnames(pNFKB_ratio_8h) = colnames(pNFAT_ratio_8h) = colnames(DP_ratio_8h) = colnames(mCherry_8h) =c("0","5","25", "50","75", "125")
######  effect size
eft <- mCherry_8h / DP_ratio_8h
eft.size <- round(t(eft[6:3,1:4]),1)[2:4,2:4]
assign(paste("eft_",i, sep=""), eft.size)

# eft123 <- rbind(eft_1,eft_2,eft_3)

pdf(file =paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS9/",
                monoclone[i],"_ABorBA_to_DP_FC.pdf",sep=""), width = 25,height = 15)

pheatmap(eft.size, 
         legend_breaks=seq(0,10,2),
         display_numbers = T,cluster_cols = F,cluster_rows=F,
         main = "effect.size",border_color = "lightgray",
         show_rownames=T, show_colnames=T,
         # labels_col = "PMA(ng/mL)",
         # labels_row = "Ionomycin(uM)",
         cellwidth = 30, cellheight =30 ,fontsize=15)

dev.off()


# display.brewer.all()
#
# color = rev(brewer.pal(9, "Reds"))
# color = colorRampPalette(color)(100)
#breaks
bk <-seq(0,100,by=1)
colors <- colorRampPalette(c( "Snow","grey61","grey51","grey11"))(length(bk))

# #####
# pdf(file =paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/FigS9/",
#                 monoclone[i],"_PI_FACS_matrix_activation_ratio.pdf",sep=""), width = 25,height = 15)

pheatmap(t(pNFAT_ratio_8h[6:3,1:4]),
        #pNFAT_activated_ratio_8h/max(pNFAT_activated_ratio_8h),
         # col= colorRampPalette(c( "GhostWhite","SkyBlue2","SteelBlue","DodgerBlue4"))(100),
        col= colors,
        legend_breaks=seq(0,100,20),
        breaks=bk,
         display_numbers = F,cluster_cols = F,cluster_rows=F,
         main = "NFAT responding ratio",border_color = "lightgray",
         show_rownames=T, show_colnames=T,
         # labels_col = "PMA(ng/mL)",
         # labels_row = "Ionomycin(uM)",
         cellwidth = 30, cellheight =30 ,fontsize=15)



pheatmap(t(pNFKB_ratio_8h[6:3,1:4]),
         #pNFKB_activated_ratio_8h/max(pNFKB_activated_ratio_8h),
         # col= colorRampPalette(c( "Snow","Wheat1","Sienna1","Sienna3"))(100),
         col= colors,
         legend_breaks=seq(0,100,20),
         breaks=bk,
         display_numbers = F,cluster_cols = F,cluster_rows=F,
         main = "NFKB responding ratio",border_color = "lightgray",
         show_rownames=T, show_colnames=T,
         # labels_col = "PMA(ng/mL)",
         # labels_row = "Ionomycin(uM)",
         cellwidth = 30, cellheight =30 ,fontsize=15)


pheatmap(t(DP_ratio_8h[6:3,1:4]),
          #pAb_activated_ratio_8h/max(pAb_activated_ratio_8h),
         col= colors,
         legend_breaks=seq(0,100,20),
         breaks=bk,
         display_numbers = F,cluster_cols = F,cluster_rows=F,
         main = "NFAT-NFKB DP responding ratio",border_color = "lightgray",
         show_rownames=T, show_colnames=T,
         # labels_col = "PMA(ng/mL)",
         # labels_row = "Ionomycin(uM)",
         cellwidth = 30, cellheight =30 ,fontsize=15)

pheatmap(t(mCherry_8h[6:3,1:4]),
        #pBa_activated_ratio_8h/max(pBa_activated_ratio_8h),
         col= colors,
        legend_breaks=seq(0,100,20),
        breaks=bk,
         display_numbers = F,cluster_cols = F,cluster_rows=F,
         main = "pBA/AB-mCherry responding ratio",border_color = "lightgray",
         show_rownames=T, show_colnames=T,
         # labels_col = "PMA(ng/mL)",
         # labels_row = "Ionomycin(uM)",
        cellwidth = 30, cellheight =30 ,fontsize=15)


dev.off()




