######################################################################################################################################################
# save.image("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_AUC_2023.Rdata")
load ("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/TCR_PI_T4h_combined_AUC_2023.Rdata")
######################################################################################################################################################

cross_0_all <- data.frame()

Treatment.TCR <- c("0.5CC", "CC", "2CC", "4CC")
for(Ti in 1:4){
  Ti_data_NFAT <- subset(All_data_NFAT_ID, Treatment.TCR == Treatment.TCR.TCR[Ti])
  Ti_data_NFKB <- subset(All_data_NFKB_ID, Treatment.TCR == Treatment.TCR.TCR[Ti])
  N <- nrow(Ti_data_NFAT)
  timelag=45
  
  AB <- list("data_NFAT"= as.data.frame(Ti_data_NFAT[,1:50]), "data_NFKB"= as.data.frame(Ti_data_NFKB[,1:50]), 
             "autocorrelation of NFAT"= data.frame(matrix(NA,N,timelag+1)),
             "autocorrelation of NFKB"= data.frame(matrix(NA,N,timelag+1)),
             "correlation of NFAT&NFKB"= data.frame(matrix(NA,N,2*timelag+1)))
  
  colnames(AB[["autocorrelation of NFAT"]])= c(0:timelag)
  colnames(AB[["autocorrelation of NFKB"]])= c(0:timelag)
  colnames(AB[["correlation of NFAT&NFKB"]])= c(-timelag:timelag)
  
  for(i2 in 1:N)
  {
    nuclNFAT_i2=as.numeric(AB[["data_NFAT"]][i2,])
    nuclNFKB_i2=as.numeric(AB[["data_NFKB"]][i2,])
    
    AB[["correlation of NFAT&NFKB"]][i2,]=ccf(nuclNFAT_i2,nuclNFKB_i2,lag.max=timelag,type = c("correlation"),plot=F)$acf
  }
  
  cross <- AB[["correlation of NFAT&NFKB"]][,paste(-timelag:timelag)]
  
  #######
  cross_0 <- data.frame(cross_correlation = cross$`0`, Treatment.TCR = rep(Treatment.TCR[Ti], nrow(cross)))
  #######
  cross_0_all <- rbind(cross_0_all, cross_0)
  #######
  
  cross_mean = apply(cross,2,mean)
  cross_errorbar=apply(cross,2,sd)
  
  assign( paste("cross_", Treatment.TCR[Ti], sep=""), cross)
  assign(paste("cross_mean_", Treatment.TCR[Ti], sep=""), cross_mean)
  assign(paste("cross_errorbar_", Treatment.TCR.TCR[Ti], sep=""), cross_errorbar)
}


library(ggplot2) 

for(Ti in 1:4){
  
  cross <- get( paste("cross_", Treatment.TCR[Ti], sep=""))
  cross_mean <- get(paste("cross_mean_", Treatment.TCR[Ti], sep=""))
  cross_errorbar <- get(paste("cross_errorbar_", Treatment.TCR[Ti], sep=""))
  
  data_plot=data.frame(auto_time=5*(as.numeric(colnames(cross))),
                       cross_mean=cross_mean,
                       cross_errorbar=cross_errorbar)
  
  pi <-  ggplot(data_plot, aes(x=auto_time, y=cross_mean)) +
    geom_line(size=1,color=rgb(126,49,142,maxColorValue=255)) +
    geom_ribbon(aes(x=auto_time,ymin=cross_mean-cross_errorbar,ymax=cross_mean+cross_errorbar),
                fill=alpha(rgb(126,49,142,maxColorValue=255),I(7/10)))+
    labs(x = "Time(min)", y = "Correlation",title = "Cross-correlation")+
    xlim(-200,200)+
    ylim(-1,1)+
    themo_demo
  
  assign(paste("p.",Treatment.TCR[Ti], sep=""), pi) 
}

# pdf(file =paste("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/z. scripts ordered by Figures/Fig1-2/Cross-correlation_2023.pdf",sep=""),width = 18,height = 8)
p.0.5CC
p.CC
p.2CC
p.4CC
# dev.off()

