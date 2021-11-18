#! /usr/bin/Rscript
library(tidyverse)
snv_path = commandArgs(trailingOnly=TRUE)[1]
ind_path = commandArgs(trailingOnly=TRUE)[2]
output   = commandArgs(trailingOnly=TRUE)[3]

snv_person.stable.id_chrom_pos_ref_alt_pass = read_tsv(snv_path,col_types=cols(.default="c")) %>% select(c(person_stable_id,chrom,pos,ref,alt,pass))
ind_person.stable.id_chrom_pos_ref_alt_pass = read_tsv(ind_path,col_types=cols(.default="c")) %>% select(c(person_stable_id,chrom,pos,ref,alt,pass))
rbind(snv_person.stable.id_chrom_pos_ref_alt_pass,ind_person.stable.id_chrom_pos_ref_alt_pass) %>% write.table(output,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

