#! /usr/bin/Rscript
# 1st arg: file 
# 2nd arg: output file name

library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(argparser)

arg_parser("funuuu") %>%
    add_argument("--denovoqc",help="funuu",default="denovoqc.20191129.txt") %>%
    add_argument("--ped",help="funuu",default="trio.invcf.sex.ped") %>%
    add_argument("--sampleqc",help="funuu",default="../sampleqc.20191202.txt") %>%
    add_argument("--snv_vq",help="funuu",default=-5.83) %>%
    add_argument("--snv_td",help="funuu",default=-5.2) %>%
    add_argument("--snv_dn",help="funuu",default=0.23) %>%
    add_argument("--snv_dg",help="funuu",default=0.02) %>%
    add_argument("--snv_df",help="funuu",flag=TRUE) %>%
    add_argument("--indel_vq",help="funuu",default=-2.86) %>%
    add_argument("--indel_td",help="funuu",default=5.5) %>%
    add_argument("--indel_df",help="funuu",flag=TRUE) %>%
    add_argument("--filteredtotal_max",help="funuu",default=11) %>%
    add_argument("--out_png",help="funuu",default="tmp.png") %>%
    add_argument("--out_txt",help="funuu",default="sampleqc.denovo.20191202.txt") %>%
    parse_args() -> argv

DENOVOQC            = argv$denovoqc
PED                 = argv$ped    
SAMPLEQC            = argv$sampleqc
SNV_VQ              = argv$snv_vq
SNV_TD              = argv$snv_td
SNV_DN              = argv$snv_dn
SNV_DG              = argv$snv_dg
SNV_DF              = if(argv$snv_df=="FALSE"){c("TRUE","FALSE") -> SNV_DF}else{"TRUE" -> SNV_DF}
INDEL_VQ            = argv$indel_vq
INDEL_TD            = argv$indel_td
INDEL_DF            = if(argv$indel_df=="FALSE"){c("TRUE","FALSE") -> INDEL_DF}else{"TRUE" -> INDEL_DF}
FILTEREDTOTAL_MAX   = argv$filteredtotal_max
OUT_PNG             = argv$out_png
OUT_TXT             = argv$out_txt

dt = read_tsv(DENOVOQC,col_types=cols(chr="c")) %>% mutate(
    vqslod          =if_else(is.na(vqslod),         min(vqslod,na.rm=T),        vqslod),
    triodenovo      =if_else(is.na(triodenovo),     min(triodenovo,na.rm=T),    triodenovo),
    dnmfilter       =if_else(is.na(dnmfilter),      min(dnmfilter,na.rm=T),     dnmfilter),
    denovogear      =if_else(is.na(denovogear),     min(denovogear,na.rm=T),    denovogear)
)
sample = read_tsv(PED,col_names=c("fid","sample","fa","mo","sex","status")) %>% filter(status==2) %>% select(sample)
qced.sample = read_tsv(SAMPLEQC) %>% mutate(trio.qc.by.snv=ifelse(
    coerced.status=="ok" & 
    contami.status=="ok" & 
    sex.sample.status!="outlier" & 
    sex.parent1.parent2=="ok" & 
    titvratio.status=="ok" & 
    nsnp.status=="ok" & 
    nins.status=="ok" & 
    ndel.status=="ok" & 
    insdelratio.status=="ok" & 
    hethomratio.status=="ok" & 
    ibd.status.parent1=="paternal" & 
    ibd.status.parent2=="paternal","ok","outlier")) %>% 
    select(sample,trio.qc.by.snv) # NA no column ga aruto kekka mo NA ninaru youda
rbind(mutate(sample,type="snv"),mutate(sample,type="indel")) -> sample_type

dt %>% group_by(sample,type) %>% summarise(count=n()) %>% right_join(sample_type) %>% 
    mutate(count=if_else(is.na(count),as.integer(0),count)) %>% ungroup %>% spread(key=type,value=count) %>% 
    mutate(total=snv+indel) %>% gather(key=type,value=count,-sample) %>% left_join(qced.sample) %>%  nest(-type) %>% 
    mutate(figure1=map(data,~ggplot(.x,aes(x=count,fill=trio.qc.by.snv)) + geom_histogram(position="stack",alpha=0.8,binwidth=1))) %>%
    mutate(figure2=map(data,~ggplot(.x,aes(x=count,fill=trio.qc.by.snv)) + geom_histogram(position="stack",alpha=0.8,binwidth=1) + xlim(-0.5,10))) -> type_data_figure

dt %>% mutate(software=case_when(
        vqslod>SNV_VQ & triodenovo>SNV_TD & dnmfilter>SNV_DN & denovogear>SNV_DG & denovofilter %in% SNV_DF & type=="snv" ~ "PASS",
        vqslod>INDEL_VQ & triodenovo>INDEL_TD & denovofilter %in% INDEL_DF & type=="indel" ~ "PASS",
        TRUE ~ "FILTERED")) %>% filter(software=="PASS") %>%
    group_by(sample,type) %>% summarise(count=n()) %>% right_join(sample_type) %>% 
    mutate(count=if_else(is.na(count),as.integer(0),count)) %>% ungroup %>% spread(key=type,value=count) %>% 
    mutate(total=snv+indel) %>% gather(key=type,value=count,-sample) %>% left_join(qced.sample) %>%  nest(-type) %>% 
    mutate(figure1=map(data,~ggplot(.x,aes(x=count,fill=trio.qc.by.snv)) + geom_histogram(position="stack",alpha=0.8,binwidth=1))) %>%
    mutate(figure2=map(data,~ggplot(.x,aes(x=count,fill=trio.qc.by.snv)) + geom_histogram(position="stack",alpha=0.8,binwidth=1) + xlim(-0.5,10))) -> filtered.type_data_figure

ggplot(dt,aes(x=type,y=vqslod,      color=sanger)) + geom_quasirandom(dodge.width=0.7,cex=2) -> g1
ggplot(dt,aes(x=type,y=triodenovo,  color=sanger)) + geom_quasirandom(dodge.width=0.7,cex=2) -> g2 
ggplot(dt,aes(x=type,y=dnmfilter,   color=sanger)) + geom_quasirandom(dodge.width=0.7,cex=2) -> g3
ggplot(dt,aes(x=type,y=denovogear,  color=sanger)) + geom_quasirandom(dodge.width=0.7,cex=2) -> g4
dt %>% group_by(type,sanger,denovofilter) %>% summarise(count=n()) %>% ggplot(aes(x=sanger,y=count,fill=denovofilter)) + facet_wrap(~type) + geom_bar(stat="identity") -> g5 #+ scale_fill_nejm() -> g5
dt %>% group_by(type,sanger,denovofilter) %>% summarise(count=n()) %>% ggplot(aes(x=sanger,y=count,fill=denovofilter)) + facet_wrap(~type) + geom_bar(stat="identity",position="fill") -> g6 #+ scale_fill_nejm() -> g6
filter(type_data_figure,type=="snv")            -> tmp; tmp$figure1[[1]] -> g7
filter(type_data_figure,type=="indel")          -> tmp; tmp$figure1[[1]] -> g8
filter(type_data_figure,type=="total")          -> tmp; tmp$figure1[[1]] -> g9
filter(type_data_figure,type=="snv")            -> tmp; tmp$figure2[[1]] -> g10
filter(type_data_figure,type=="indel")          -> tmp; tmp$figure2[[1]] -> g11
filter(type_data_figure,type=="total")          -> tmp; tmp$figure2[[1]] -> g12
filter(filtered.type_data_figure,type=="snv")   -> tmp; tmp$figure1[[1]] -> g13
filter(filtered.type_data_figure,type=="indel") -> tmp; tmp$figure1[[1]] -> g14
filter(filtered.type_data_figure,type=="total") -> tmp; tmp$figure1[[1]] -> g15
filter(filtered.type_data_figure,type=="snv")   -> tmp; tmp$figure2[[1]] -> g16
filter(filtered.type_data_figure,type=="indel") -> tmp; tmp$figure2[[1]] -> g17
filter(filtered.type_data_figure,type=="total") -> tmp; tmp$figure2[[1]] -> g18

plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,nrow=3)
ggsave(OUT_PNG,width=50,height=20,limitsize = FALSE)

inner_join(
filtered.type_data_figure %>% unnest(data) %>% spread(key=type,value=count) %>% rename(snv.filtered=snv,indel.filtered=indel,total.filtered=total),
type_data_figure %>% unnest(data) %>% spread(key=type,value=count)) %>% 
mutate(total.filtered.sampleqc=if_else(total.filtered>=FILTEREDTOTAL_MAX,"outlier","ok")) %>% write_tsv(OUT_TXT)
