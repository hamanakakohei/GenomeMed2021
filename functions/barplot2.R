#! /usr/bin/Rscript
# 1st arg: 1st col: sample 2nd col: depth
# 2nd arg: output png file name

library(tidyverse)
FILE = commandArgs(trailingOnly=TRUE)[1]
OUTPUT = commandArgs(trailingOnly=TRUE)[2]

PC = read_tsv(FILE)
png(OUTPUT)
PC$Total_Depth %>% as.vector %>% barplot
dev.off()
