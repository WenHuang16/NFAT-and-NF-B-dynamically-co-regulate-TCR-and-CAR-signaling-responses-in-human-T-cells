library (VennDiagram)

#####
input_dir <- "/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/CAR/rawcounts"

setwd(input_dir)

ExpreFile_basal <- dir(paste(input_dir,"basal_raw",sep="/"))
ExpreFile_activation <- dir(paste(input_dir,"activation_raw",sep="/"))
ExpreFile_act_NC <- dir(paste(input_dir,"act_NC",sep="/"))
type <-  c("basal", "activation","act_NC")          

###################################################################################################################
library(DESeq2)

for (typ in 1:3){
  sample <- get(paste("ExpreFile", type[typ],sep="_"))
  
  sample.frame <- data.frame(matrix(0,length(sample),5))
  colnames(sample.frame) <- c("samplename","treatment", "timepoint", "rep","group.name")
  
  for (x in 1:length(sample))
  {
    a <- strsplit(sample[x],split="[-/_]")
    sample.frame[x,1] <- paste(a[[1]][1],a[[1]][2],a[[1]][3],sep="-")
    sample.frame[x,2] <- a[[1]][1]
    sample.frame[x,3] <- a[[1]][2]
    sample.frame[x,4] <- a[[1]][3]
    sample.frame[x,5] <- paste(a[[1]][1],a[[1]][2],sep=".")
  }
  assign(paste("sample.frame", type[typ],sep="_"), sample.frame)
  
  
  sampleTable.typ	<- data.frame(samplename = sample.frame$samplename,
                              filename = sample,
                              condition = factor(sample.frame$group.name),
                              cartype = sample.frame$treatment)
  assign(paste("sampleTable",type[typ],sep="_"), sampleTable.typ)
}

# ############################################# CAR basal #############################################
# 
#   dds <- DESeqDataSetFromHTSeqCount(sampleTable_basal, directory = "./basal_raw", design= ~ condition)
#   dds$condition <- relevel(dds$condition, ref = "J.NC")
#   
#   dds2 <- DESeq(dds)
# 
# 
# cont <- resultsNames(dds2)
# 
# tres <- c("J96.NC","J97.NC","J99.NC","J104.NC")
# DEG.all_basal.up <- data.frame(matrix(2,1,0))
# 
# 
# for( tre in 1:4){
#     name <- paste("condition_",tres[tre],"_vs_J.NC",sep="")
#     res <- results(dds2, name = name)
# 
#     # sort by p.adj
#     res <- res[order(res$pvalue),]
#     # get diff_gene
#     diff_gene <- subset(as.data.frame(res), pvalue < 0.05 & log2FoldChange > 0)
# 
#     ################# symbol conversion
#    library(org.Hs.eg.db)
#     idfound <- rownames(diff_gene) %in% mappedRkeys(org.Hs.egREFSEQ)
#     exprSet_annotated <- data.frame(diff_gene[idfound,])
# 
#     REFSEQ <- toTable(org.Hs.egREFSEQ)
#     SYMBOL <- toTable(org.Hs.egSYMBOL)
# 
#     m <- match(rownames(exprSet_annotated), REFSEQ$accession)
#     exprSet_annotated$EntrezGene <- REFSEQ$gene_id[m]
# 
#     m2 <- match(exprSet_annotated$EntrezGene, SYMBOL$gene_id)
#     exprSet_annotated$symbol <- SYMBOL$symbol[m2]
#     
#     DEG.all_basal.up <-  unique(rbind(DEG.all_basal.up, data.frame(RefID = rownames(exprSet_annotated), symbol=exprSet_annotated$symbol)))
# 
# }
# 
# ######## get FC list ########
# list.DEG.all_basal.up <- list()
# for( tre in 1:4){
#   name <- paste("condition_",tres[tre],"_vs_J.NC",sep="")
#   res <- results(dds2, name = name)
# 
#   # get diff_gene in DEG.syb.all_basal.up
#   x = as.data.frame(res)
#   diff_gene <- subset(x, rownames(x) %in% DEG.all_basal.up$RefID)
#   diff_gene$RefID <- rownames(diff_gene)
# 
#   list.DEG.all_basal.up[[tre]] <- merge(diff_gene,DEG.all_basal.up, by="RefID")
# }
#   names(list.DEG.all_basal.up) <- tres
# 
# 
# ############################################# CAR ACTIVATION #############################################
# 
# DEG.all_act.up <- data.frame(matrix(2,1,0))
#   
# for (ci in 1:4){
#   
#   car.type <- c("J96","J97","J99","J104")
#   
#   sampleTable.ci	<- subset(sampleTable_activation, sampleTable_activation$cartype == car.type[ci])
#   
#   dds <- DESeqDataSetFromHTSeqCount(sampleTable.ci, directory = "./activation_raw", design= ~ condition)
#   dds$condition <- relevel(dds$condition, ref = paste(car.type[ci], "NC", sep="."))
#   dds2 <- DESeq(dds)
#   
#   # sort by p.adj
#   res <- results(dds2)    
#   res <- res[order(res$pvalue),]
#       
#       # get diff_gene
#       diff_gene <- subset(as.data.frame(res), pvalue < 0.05 & log2FoldChange > 0)
#       
#       ################# symbol conversion
#       library(org.Hs.eg.db)
#       idfound <- rownames(diff_gene) %in% mappedRkeys(org.Hs.egREFSEQ)
#       exprSet_annotated <- data.frame(diff_gene[idfound,])
#       
#       REFSEQ <- toTable(org.Hs.egREFSEQ)
#       SYMBOL <- toTable(org.Hs.egSYMBOL)
#       
#       m <- match(rownames(exprSet_annotated), REFSEQ$accession)
#       exprSet_annotated$EntrezGene <- REFSEQ$gene_id[m]
#       
#       m2 <- match(exprSet_annotated$EntrezGene, SYMBOL$gene_id)
#       exprSet_annotated$symbol <- SYMBOL$symbol[m2]
# 
#       DEG.all_act.up <-  unique(rbind(DEG.all_act.up, data.frame(RefID = rownames(exprSet_annotated), symbol=exprSet_annotated$symbol)))
#       
# }
#   
# ######## get CAR activation FC list ########
# list.DEG.all_act.up <- list()
# for( ci in 1:4){
#   car.type <- c("J96","J97","J99","J104")
#   sampleTable.ci	<- subset(sampleTable_activation, sampleTable_activation$cartype == car.type[ci])
#   
#   dds <- DESeqDataSetFromHTSeqCount(sampleTable.ci, directory = "./activation_raw", design= ~ condition)
#   dds$condition <- relevel(dds$condition, ref = paste(car.type[ci], "NC", sep="."))
#   dds2 <- DESeq(dds)
#   
#   # sort by p.adj
#   res <- results(dds2)    
#     
#   # get diff_gene in DEG.syb.all_basal.up
#   x = as.data.frame(res)
#   diff_gene <- subset(x, rownames(x) %in% DEG.all_act.up$RefID)
#   diff_gene$RefID <- rownames(diff_gene)
#     
#   list.DEG.all_act.up[[ci]] <- merge(diff_gene,DEG.all_act.up, by="RefID")
# }
# names(list.DEG.all_act.up) <- car.type
#   

############################################# JSH VS ACTIVATED CAR #############################################
############################################# TO COMPARE TCR ACTIVATION AND CAR ACTIVATION #############################################

dds <- DESeqDataSetFromHTSeqCount(sampleTable_act_NC, directory = "./act_NC", design= ~ condition)
dds$condition <- relevel(dds$condition, ref = "J.NC")

dds2 <- DESeq(dds)
cont <- resultsNames(dds2)

tres <- c("J96.treat","J97.treat","J99.treat","J104.treat","J.treat")
DEG.all_actCAR.up <- data.frame(matrix(2,1,0))

for( tre in 1:5){
  # for (tp in 1:3){
    name <- paste("condition_",tres[tre],"_vs_J.NC",sep="")
    res <- results(dds2, name = name)
    
    # sort by p.adj
    res <- res[order(res$pvalue),]
    
    # get diff_gene
    diff_gene <- subset(as.data.frame(res), pvalue < 0.05 & log2FoldChange > 0)
    
    ################# symbol conversion
    library(org.Hs.eg.db)
    idfound <- rownames(diff_gene) %in% mappedRkeys(org.Hs.egREFSEQ)
    exprSet_annotated <- data.frame(diff_gene[idfound,])
    
    REFSEQ <- toTable(org.Hs.egREFSEQ)
    SYMBOL <- toTable(org.Hs.egSYMBOL)
    
    m <- match(rownames(exprSet_annotated), REFSEQ$accession)
    exprSet_annotated$EntrezGene <- REFSEQ$gene_id[m]
    
    m2 <- match(exprSet_annotated$EntrezGene, SYMBOL$gene_id)
    exprSet_annotated$symbol <- SYMBOL$symbol[m2]
    
    # resdata <- merge(as.data.frame(exprSet_annotated),
    #                  as.data.frame(log2(counts(dds, normalized=F))),
    #                  by="row.names", sort=FALSE)
    # 
    # assign(paste("DEG.actNC.",tres[tre],sep=""),resdata)
    
    DEG.all_actCAR.up <- unique(rbind(DEG.all_actCAR.up, data.frame(RefID = rownames(exprSet_annotated), symbol=exprSet_annotated$symbol)))
    
  # }
}

#################### get CAR/TCR activation FC list ########
list.DEG.all_actCAR.up <- list()

for( tre in 1:5){
  name <- paste("condition_",tres[tre],"_vs_J.NC",sep="")
  res <- results(dds2, name = name)
  
  # get diff_gene in DEG.all_actCAR.up
  x = as.data.frame(res)
  diff_gene <- subset(x, rownames(x) %in% DEG.all_actCAR.up$RefID)
  diff_gene$RefID <- rownames(diff_gene)
  
  list.DEG.all_actCAR.up[[tre]] <- merge(diff_gene, DEG.all_actCAR.up, by="RefID")
}
names(list.DEG.all_actCAR.up) <- tres

###############################
#  DEG.all_actCAR_dataframe  ##
###############################

DEG.all_actCAR_dataframe <- data.frame(symbol= list.DEG.all_actCAR.up[[1]]$symbol, RefID=list.DEG.all_actCAR.up[[1]]$RefID)

for( tre in 1:5){ ###  
  DEG.all_actCAR_dataframe[,tre+2] <- list.DEG.all_actCAR.up[[tre]]$log2FoldChange
}
colnames(DEG.all_actCAR_dataframe) <- c("symbol","RefID",paste("log2FC",names(list.DEG.all_actCAR.up),sep="."))

save.image("/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/Fig4_RNA_DEseq2.by.group_CAR_2023.RData") # produced by Fig4_RNA_DESeq2.by.group.CAR




