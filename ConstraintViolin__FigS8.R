library(tidyverse)

gene_new_enst__path		="gene_new_enst.txt"
clingen_path	   		="ClinGen_HI.txt"
MisEnr26_path	   		="MisEnr26Gene.txt"
ENST_OF_INTEREST__path		="enst.coding.canonical.gencode.noalthaplo.noenstr.txt"
loeuf_path         		="gnomad.v2.1.1.lof_metrics.by_transcript.txt"
gene_enst__ddg2p.ad.hi__path	="gene_enst__ddg2p.ad.hi.txt"

gene__NotLof=c("PSMC3","CBX5","PIP5K1C","MAST3","SUPT16H","SEPT2","KBTBD7","HIST1H2AE","BRD3")
enst__NotLof=c("ENST00000379483","ENST00000303910","ENST00000216297","ENST00000303407","ENST00000589578","ENST00000298852","ENST00000401990","ENST00000209875","ENST00000262811")
read_tsv(ENST_OF_INTEREST__path,col_names="enst")                                                                                                                                            -> enst__of.interest
read_tsv(gene_enst__ddg2p.ad.hi__path)                                                                                                                                                      -> gene_enst__ddg2p.ad.hi
read_tsv(clingen_path,col_names="gene") %>% inner_join(enst_gene) %>% inner_join(enst__of.interest) %>% select(enst) %>% mutate(ClinGenHI="yes") -> enst_ClinGenHI
read_tsv(loeuf_path) %>% dplyr::select(c(transcript,pLI,oe_lof,oe_lof_upper,oe_mis,mis_z)) %>% rename(enst=transcript,oe.mis=oe_mis,mis.z=mis_z) -> enst_pLI_oe_lof_oe_lof_upper
read_tsv(gene_new_enst__path,col_types="cdc",col_names=c("gene","new","enst")) %>% mutate(new=ifelse(new==1,"new","old"),new=ifelse(enst %in% enst__NotLof,"NotLof",new))-> new_enst
gene_enst__ddg2p.ad.hi %>% mutate(ddg2p="yes") %>% select(-gene) -> enst_ddg2p
enst__of.interest %>% 
    left_join(enst_pLI_oe_lof_oe_lof_upper) %>% 
    left_join(new_enst) %>% 
    left_join(enst_ddg2p) %>% 
    left_join(enst_ClinGenHI) %>% 
    filter(!is.na(pLI) & !is.na(oe_lof) & !is.na(oe_lof_upper)) %>% 
    distinct(enst,.keep_all=TRUE) %>% 
    mutate(new=ifelse(is.na(new) & is.na(ddg2p) & is.na(ClinGenHI),"control",new)) %>% 
    filter(new %in% c("new","control")) %>% 
    ggplot(aes(x=new,y=pLI,fill=new)) + geom_violin()

read_tsv(MisEnr26_path,col_names="gene") %>% inner_join(new_enst) %>% pull(enst) -> mis.enst__v
read_tsv(gene_new_enst__path,col_types="cdc",col_names=c("gene","new","enst")) %>% mutate(new=ifelse(new==1,"new","old"),new=ifelse(enst %in% mis.enst__v,new,"NotMis"))-> new_enst
enst__of.interest %>% 
    left_join(enst_pLI_oe_lof_oe_lof_upper) %>% 
    left_join(new_enst) %>% 
    left_join(enst_ddg2p) %>% 
    left_join(enst_ClinGenHI) %>% 
    filter(!is.na(pLI) & !is.na(oe_lof) & !is.na(oe_lof_upper) & !is.na(oe.mis) & !is.na(mis.z)) %>% 
    distinct(enst,.keep_all=TRUE) %>% 
    mutate(new=ifelse(is.na(new) & is.na(ddg2p) & is.na(ClinGenHI),"control",new)) %>% 
    filter(new %in% c("new","control")) %>% 
    ggplot(aes(x=new,y=oe.mis,fill=new)) + geom_violin() + ylim(-8,NA) #+ ylim(-2,NA)

