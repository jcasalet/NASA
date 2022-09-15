library(DESeq2)
library("vsn")

args <- commandArgs(trailingOnly = TRUE)
print(args)
exprFile <- args[1]
metaFile <- args[2]
morFile <- args[3]

#exprFile <- '/Users/jcasalet/Desktop/RESEARCH/LIVER/HETEROSKEDASCTICITY/S2N/FINAL/R_DESEQ2/expr_unnorm_protein-coding_top-0_gene_x_sample.csv'
#metaFile <- '/Users/jcasalet/Desktop/RESEARCH/LIVER/HETEROSKEDASCTICITY/S2N/FINAL/R_DESEQ2/metadata.csv'
#morFile <- '/Users/jcasalet/Desktop/RESEARCH/LIVER/HETEROSKEDASCTICITY/S2N/FINAL/R_DESEQ2/mor.csv'
#batch <- 'dataset'

expr <- read.csv(exprFile, header=TRUE, row.names=1, stringsAsFactors=TRUE, check.names=FALSE)
meta <- read.csv(metaFile, header=TRUE, row.names=1, stringsAsFactors=TRUE, check.names=FALSE)
#colData <- DataFrame(condition=factor(meta$oro_thresh), libprep=factor(meta$dataset))
colData <- DataFrame(condition=factor(meta$oro_thresh), libprep=factor(meta$libprep))
ncol(expr)
nrow(colData)

if (nlevels(colData$libprep) != 1) {
  dds <- DESeqDataSetFromMatrix(ceiling(expr), colData=colData, ~libprep)
} else {
  dds <- DESeqDataSetFromMatrix(ceiling(expr), colData=colData, ~1)
}

dds <- DESeq(dds)
dds_sizefactors <- estimateSizeFactors(dds)
dds_normalized <- counts(dds_sizefactors, normalize=TRUE)
#write.csv(as.data.frame(dds_normalized, check.names=FALSE), file=morFile, quote=FALSE)
#write.csv(data.frame(dds_normalized, check.names=FALSE), file=morFile, quote=FALSE)
write.csv(data.frame("gene"=rownames(dds_normalized), dds_normalized, check.names=FALSE), file=morFile, quote=FALSE)
