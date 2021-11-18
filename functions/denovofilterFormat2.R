#! /usr/bin/Rscript
library(tidyverse)
varid_pt_fa_mo_path     = commandArgs(trailingOnly=TRUE)[1]
varid_sample_depth_path = commandArgs(trailingOnly=TRUE)[2]
id_consequence_symbol_enst       = commandArgs(trailingOnly=TRUE)[3]
outputsnv               = commandArgs(trailingOnly=TRUE)[4]
outputindel             = commandArgs(trailingOnly=TRUE)[5]

varid_pt_fa_mo = read_tsv(varid_pt_fa_mo_path,col_names=c("varid","pt","fa","mo"))
varid_sample_depth = read_delim(varid_sample_depth_path,col_names=c("varid","sample","depth"),delim=" ")
chrom_pos_ref_alt_consequence_symbol = read_tsv(id_consequence_symbol_enst,col_names=c("id","consequence","symbol")) %>% separate(id,c("chrom","pos","ref","alt"),"-")

varid_pt_fa_mo %>% select(-varid) %>% mutate(fam = pt) %>% distinct(fam,.keep_all=TRUE) %>% 
    gather(key=relation,value=sample,pt,fa,mo) -> fam_relation_sample

varid_sample_depth %>% left_join(fam_relation_sample) %>% select(-sample) %>% spread(key=relation,value=depth) %>%
    separate("varid",into=c("chrom","pos","ref","alt"),sep="-") %>%
    rename(dp4_child=pt,dp4_father=fa,dp4_mother=mo,person_stable_id=fam) %>% 
    mutate(max_af=0,type=if_else(str_length(ref)==1 & str_length(alt)==1,"snv","indel")) %>%
    left_join(chrom_pos_ref_alt_consequence_symbol,by=c("chrom","pos","ref","alt")) -> dt

filter(dt,type=="indel") %>% select(-type) %>% select(c(person_stable_id,chrom,pos,ref,alt,symbol,consequence,dp4_child,dp4_father,dp4_mother,max_af)) %>% write_tsv(outputindel)
filter(dt,type=="snv") %>% select(-type) %>% mutate(pp_dnm=1,in_child_vcf=1,in_mother_vcf=0,in_father_vcf=0) %>% select(c(person_stable_id,chrom,pos,ref,alt,symbol,consequence,dp4_child,dp4_father,dp4_mother,max_af,pp_dnm,in_child_vcf,in_mother_vcf,in_father_vcf)) %>% write_tsv(outputsnv)

