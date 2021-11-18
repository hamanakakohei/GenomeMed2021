#! /usr/bin/Rscript
# 1st arg: file (header: Variant_ID Sample_ID1 Sample_ID2,,,; 1st row: Variant_ID, Variant_ID1, Vairnat_ID2,,,
# 1st arg: file (1st row: "Variant"
# 2nd arg: n of pc for plot
# 3rd arg: output plot name
# 4th arg: output table name

library(tidyverse)
library(stats)
FILE1 = commandArgs(trailingOnly=TRUE)[1]
N.PC.PLOT = commandArgs(trailingOnly=TRUE)[2]
OUTPUT.PLOT = commandArgs(trailingOnly=TRUE)[3]
OUTPUT.TABLE = commandArgs(trailingOnly=TRUE)[4]

VARIANT_SAMPLES = read_tsv(FILE1)
VARIANT_SAMPLES[,-1] %>% as.matrix %>% t %>% scale %>% prcomp -> RES

RES$x %>% write.table(OUTPUT.TABLE,quote=F,sep="\t")

png(OUTPUT.PLOT)
summary(RES)$importance[2,] %>% as.vector %>% head(n=N.PC.PLOT) %>% barplot 
dev.off()

