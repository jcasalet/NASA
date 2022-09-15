library(DESeq2)
library("vsn")

args <- commandArgs(trailingOnly = TRUE)
print(args)
exprFile <- args[1]
metaFile <- args[2]
outFile <- args[3]
transformation <- args[4]

#exprFile <- '/Users/jcasalet/Desktop/RESEARCH/LIVER/DATA/JC/SOUP2NUTS/expr_before_shrink.csv'
#metaFile <- '/Users/jcasalet/Desktop/RESEARCH/LIVER/DATA/JC/SOUP2NUTS/meta.csv'
#transformation <- "vsd"

expr <- read.csv(exprFile, header=TRUE, row.names=1, sep=',', check.names=FALSE)
meta <- read.csv(metaFile, header=TRUE, row.names=1, sep=',', check.names=FALSE)

dim(expr)
head(expr)
colData <- DataFrame(condition=factor(meta$oro_thresh), libprep=factor(meta$libprep))
contrast = c("condition", "0", "1")
dds <- DESeqDataSetFromMatrix(ceiling(expr), colData=colData, ~condition + libprep)
#keep <- rowSums(counts(dds) >= 5) > 10
#dds <- dds[keep,]

dds <- DESeq(dds)

if(transformation == 'vsd')
{
  vsd <- vst(dds, blind=TRUE)
  meanSdPlot(assay(vsd))
  vsd_df <- assay(vsd)
  write.csv(data.frame(vsd_df, check.names=FALSE), file=outFile)
} else if(transformation == 'rld')
{
  rld <- rlog(dds, blind=TRUE)
  meanSdPlot(assay(rld))
  rld_df <- assay(rld)
  write.csv(data.frame(rld_df, check.names=FALSE), file=outFile)
} else if(transformation == 'ntd')
{
  ntd <- normTransform(dds)
  meanSdPlot(assay(ntd))
  ntd_df <- assay(ntd)
  write.csv(data.frame(ntd_df, check.names=FALSE), file=outFile)
} 
