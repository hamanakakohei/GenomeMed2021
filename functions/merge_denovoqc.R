#! /usr/bin/Rscript

library(tidyverse)
FILE1 = commandArgs(trailingOnly=TRUE)[1]
FILE2 = commandArgs(trailingOnly=TRUE)[2]
FILE3 = commandArgs(trailingOnly=TRUE)[3]
FILE4 = commandArgs(trailingOnly=TRUE)[4]
FILE5 = commandArgs(trailingOnly=TRUE)[5]
FILE6 = commandArgs(trailingOnly=TRUE)[6]
FILE7 = commandArgs(trailingOnly=TRUE)[7]
FILE8 = commandArgs(trailingOnly=TRUE)[8]
FILE9 = commandArgs(trailingOnly=TRUE)[9]
OUTPUT = commandArgs(trailingOnly=TRUE)[10]

chr_pos_ref_alt_sample_fa_mo = read_tsv(FILE1,col_types=cols(.default="c"),col_names=c("id","sample","fa","mo")) %>% separate(id,c("chr","pos","ref","alt"),sep="-")
chr_pos_ref_alt_vq           = read_tsv(FILE2,col_types=cols(.default="c"),col_names=c("chr","pos","ref","alt","vq"))
sample_chr_pos_ref_alt_td    = read_tsv(FILE3,col_types=cols(.default="c"),col_names=c("sample","chr","pos","ref","alt","td"))
sample_chr_pos_dn            = read_tsv(FILE4,col_types=cols(.default="c"),col_names=c("sample","chr","pos","ref","alt","dn")) %>% select(-c(ref,alt))
sample_chr_pos_ref_alt_dg    = read_tsv(FILE5,col_types=cols(.default="c"),col_names=c("sample","chr","pos","ref","alt","dg"))
sample_chr_pos_ref_alt_df    = read_tsv(FILE6,col_types=cols(.default="c"),col_names=c("sample","chr","pos","ref","alt","df"))
sample_chr_pos_ref_alt_sa    = read_tsv(FILE7,col_types=cols(.default="c")) %>% mutate(sa="confirmed")
chr_pos_ref_alt_snpeff       = read_tsv(FILE8,col_types=cols(.default="c"),col_names=c("id","snpeff.csq","snpeff.gene","snpeff.enst","snpeff.variant")) %>% separate(id,c("chr","pos","ref","alt"),sep="-")
chr_pos_ref_alt_annovar      = read_tsv(FILE9,col_types=cols(.default="c")) %>% separate(`chr-pos-ref-alt`,c("chr","pos","ref","alt"),sep="-")

chr_pos_ref_alt_sample_fa_mo %>%
    mutate(type=if_else(str_length(ref)==1 & str_length(alt)==1,"snv","indel")) %>%
    left_join(chr_pos_ref_alt_vq) %>%
    left_join(sample_chr_pos_ref_alt_td) %>%
    left_join(sample_chr_pos_dn) %>%
    left_join(sample_chr_pos_ref_alt_dg) %>%
    left_join(sample_chr_pos_ref_alt_df) %>%
    left_join(sample_chr_pos_ref_alt_sa) %>%
    left_join(chr_pos_ref_alt_snpeff) %>%
    left_join(chr_pos_ref_alt_annovar) %>%
    rename(vqslod=vq,triodenovo=td,dnmfilter=dn,denovogear=dg,denovofilter=df,sanger=sa) %>%
    write_tsv(OUTPUT)
