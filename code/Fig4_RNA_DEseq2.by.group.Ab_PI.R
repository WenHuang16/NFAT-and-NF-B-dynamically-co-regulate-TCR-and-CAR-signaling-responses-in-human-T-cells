
#####
input_dir <- "/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/rawcounts/"
setwd(input_dir)
ExpreFile <- dir(input_dir)

###################################################################################################################

library(DESeq2)

sample <- ExpreFile
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



sampleTable	<- data.frame(samplename = sample.frame$samplename,
                          filename = ExpreFile,
                          condition = factor(sample.frame$group.name))

dds <- DESeqDataSetFromHTSeqCount(sampleTable, directory = ".", design= ~ condition)

dds$condition <- relevel(dds$condition, ref = "NC.0")

dds2 <- DESeq(dds)
cont <- resultsNames(dds2)

tres <- c("P","I","PIC","S")
tps <- c("030","120","240")
DEG.all_act.up <- data.frame()
# DEG.all.down <- data.frame()

for( tre in 1:4){ ###  consider soluble antibody here
  for (tp in 1:3){
    name <- paste("condition_",tres[tre],".",tps[tp],"_vs_NC.0",sep="")
    res <- results(dds2, name = name)

    # sort by p.adj
    res <- res[order(res$padj),]

    # get diff_gene
    diff_gene <- subset(as.data.frame(res), padj < 0.05 & log2FoldChange > 0)

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


    resdata <- merge(as.data.frame(exprSet_annotated),
                     as.data.frame(counts(dds2, normalized=T)),
                     by="row.names", sort=FALSE)

    assign(paste("DEG.",tres[tre],".",tps[tp],sep=""),resdata)

    DEG.all_act.up <-  unique(rbind(DEG.all_act.up, data.frame(RefID = rownames(exprSet_annotated), symbol=exprSet_annotated$symbol)))
    # DEG.all.down <- unique(rbind(DEG.all.down, data.frame(exprSet_annotated$symbol)))
  }
}


# DEGs <- unique(rbind(DEG.all.up, DEG.all.down))

write.csv(DEG.all.up,file = "./DEG.all.up.PorIorPI.csv",row.names = F)

############################################### extract Fold change of DEGs in each condition ###########################################
######## get CAR activation FC list ########
list.DEG.all_act.up <- list()
i=0

for( tre in 1:4){ ###  consider soluble antibody here
  for (tp in 1:3){
    i=i+1
    name <- paste("condition_",tres[tre],".",tps[tp],"_vs_NC.0",sep="")
    res <- results(dds2, name = name)

    # get diff_gene in DEG.syb.all_basal.up
    x = as.data.frame(res)
    diff_gene <- subset(x, rownames(x) %in% DEG.all_act.up$RefID)
    diff_gene$RefID <- rownames(diff_gene)


    list.DEG.all_act.up[[i]] <- merge(diff_gene,DEG.all_act.up, by="RefID")
  }
}

###########
DEG.all_act.up_dataframe <- data.frame(RefID = list.DEG.all_act.up[[1]]$RefID, symbol=list.DEG.all_act.up[[1]]$symbol)
namei <- data.frame()
i=0
for( tre in 1:4){ ###  consider soluble antibody here
  for (tp in 1:3){
    namei <- rbind(namei, paste(tres[tre],".",tps[tp],"_vs_NC.0",sep=""))

    i=i+1
    DEG.all_act.up_dataframe[,i+2] <- list.DEG.all_act.up[[i]]$log2FoldChange
  }
}

names(list.DEG.all_act.up) <- namei$X.P.030_vs_NC.0.
colnames(DEG.all_act.up_dataframe) <- c("RefID", "Symbol", namei$X.P.030_vs_NC.0.)


# save.image("/storage/labnas5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/Fig4_RNA_DEseq2.by.group_Ab_PI_2023.RData")

# load("//LINLAB5/LabNAS-5/HuangWen2023/T cell dynamics paper/Sequencing/metadata/Fig4_RNA_DEseq2.by.group_Ab_PI_2023.RData")
