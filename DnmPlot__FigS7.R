library(tidyverse)
library(zoo)
library(cowplot)
res__path                           = "res20210302withoutCNV.rds"
chr_start_end_enst__path            = "enst_start_end.txt"
depth__path                         = "gnomad.exomes.r2.1.coverage.tsv.bgz"
pfam__path                          = "pfam.ucscgene.hg19.gtf"
gene.region_path                    = "gencode.v19.annotation.gff3.gz"
variant_infos__gnomad__path         = "variant_infos__gnomad.txt.gz"
enst_class__path		    = "enst__NewCandidate.txt"
MIS                                 = c("missense_variant","missense_variant&splice_region_variant","protein_protein_contact","rare_amino_acid_variant")

read_tsv(gene.region_path,comment="#",col_types="c_cdd_c_c",col_names=c("chr","type","start","end","strand","att")) %>%
    separate(att,sep=";",c("att1","att2","att3","att4")) %>%
    separate(att4,sep="=",c("n","enst")) %>%  
    separate(enst,sep="\\.",c("enst")) %>%
    mutate(size=end-start+1) %>%
    filter(chr!="chrM" & type=="CDS") %>%
    select(-c(att1,att2,att3,n,type)) -> chr_start_end_strand_enst_size

read_tsv(enst_class__path) -> enst_class__MyInterest
read_tsv(variant_infos__gnomad__path,col_types=cols(chr="c")) %>% separate(ANN,c("n1","csq","n2","gene","n3","n4","enst","n5","n6","n7","aa"),sep="\\|") %>%
    separate(enst,c("enst","n5"),":") %>% filter(filter=="PASS" & csq %in% MIS) %>% select(c(chr,pos,gene,enst,aa)) %>% mutate(chr=paste0("chr",chr))                       -> chr_pos_gene_enst_aa__missense.gnomad
read_tsv(pfam__path,col_types="c__dd___c",col_names=c("chr","start","end","pfam")) %>% separate(pfam,c("n","pfam")," ") %>% separate(pfam,c("n","pfam"),'"') %>% select(-n) -> chr_start_end_pfam
read_tsv(gzfile(depth__path),col_types="_c_d") %>% separate(locus,c("chr","pos")) %>% mutate(chr=paste0("chr",chr),pos=as.numeric(pos)) %>% rename(depth.gnomad=median)     -> chr_pos_depth.gnomad
read_tsv(chr_start_end_enst__path,col_names=c("chr","start","end","enst"))                                                                                                  -> chr_start_end_enst
readRDS(res__path)                                                                                                                                                          -> enrichment.res
chr_start_end_pfam %>% mutate(df__all.base=map2(start,end,~seq(.x,.y,1) %>% enframe(value="pos") %>% select(-name))) %>% unnest(df__all.base) %>% select(c(chr,pos,pfam))   -> chr_pos_pfam
chr_pos_gene_enst_aa__missense.gnomad %>% group_by(chr,pos,enst) %>% summarise(count.aa=n())                                                                                -> chr_pos_enst_count.aa
enrichment.res %>% inner_join(enst_class__MyInterest) %>% mutate(map2(enst,class,plot_fig))


# functions
gather_overlap_domains = function(DF){ nest(DF,pfam) %>% mutate(data=map(data,~.x %>% arrange %>% pull(pfam) %>% paste(collapse=";"))) %>% unnest %>% rename(pfam=data) %>% mutate(pfam=ifelse(pfam=="NA",NA,pfam)) %>% return }

plot_fig = function(enst__of.interest, enrichment.type){
    OUT = paste0(enst__of.interest,"-",enrichment.type,".png")
    chr_start_end_strand_enst_size %>% filter(enst==enst__of.interest) %>% 
        mutate(df__all.base=map2(start,end,~seq(.x,.y,1) %>% enframe(value="pos") %>% select(-name))) %>%
        unnest(df__all.base) %>% add_pos.aa %>% add_n.exon %>% left_join(chr_pos_depth.gnomad) %>% left_join(chr_pos_pfam) %>% gather_overlap_domains %>% left_join(chr_pos_enst_count.aa) %>% 
        mutate(count.aa=ifelse(is.na(count.aa),0,count.aa)) %>% mutate(count.aa.rollmean=rollmean(count.aa,15,fill=0)) -> data
    
    do.call("rbind", enrichment.res %>% filter(enst==enst__of.interest & class==enrichment.type) %>% pull(meta) %>% .[[1]] ) %>% extract_sample2 %>% filter(!is.na(anno)) -> anno_pos_id
    plot_fig2(enst__of.interest,data,anno_pos_id,OUT)
}

plot_fig2 = function(ENST,DF,anno_pos_id,OUT){
    g__theme = theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), plot.margin=unit(c(0,0,0,0),"cm"))
    DF$pos.aa %>% max -> POS.AA.MAX
    DF %>% group_by(n.exon) %>% summarise(pos.aa.start=min(pos.aa), pos.aa.end=max(pos.aa)) %>% mutate(y__min=0,y__max=1) -> DF__for.exon
    DF %>% group_by(pos.aa) %>% summarise(count.aa.rollmean=mean(count.aa.rollmean), depth.gnomad=mean(depth.gnomad)) -> DF2
    DF2 %>% ggplot(aes(x=pos.aa, y=count.aa.rollmean, color="gray")) + geom_line() + geom_area(fill="gold3",color="gold3") + g + g__theme + ylim(0,0.5) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())  -> g1
    DF2 %>% ggplot(aes(x=pos.aa, y=depth.gnomad     , color="gray")) + geom_line() + geom_area(fill="lightblue",color="lightblue") + g + g__theme + ylim(0,100) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())   -> g2
    DF__for.exon %>% ggplot() + geom_rect(aes(xmin=pos.aa.start,xmax=pos.aa.end,ymin=y__min,ymax=y__max), color="black", fill="gray85") + g + g__theme + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())-> g3
    DF %>% mutate(pfam=ifelse(is.na(pfam),"nodomain",pfam)) %>% annotate_pfam %>% ggplot() + geom_rect(aes(xmin=pos.aa.start,xmax=pos.aa.end,ymin=y__min,ymax=y__max,fill=pfam)) + g + g__theme + theme(legend.position="none") + xlim(0, POS.AA.MAX) -> g4
    anno_pos_id %>% mutate(
        study=case_when(
            str_detect(id, "(DDD|GDX|RUMC)") ~ "DDD",
            str_detect(id, "YCU") ~ "YCU",
            TRUE ~ "denovodb"
        ), 
        HEIGHT=case_when(
            study == "DDD"      ~ 0,
            study == "YCU"      ~ -0.06,
            study == "denovodb" ~ -0.03
        ),
        COLOR=case_when(
            study == "DDD"      ~ "red",
            study == "YCU"      ~ "blue",
            study == "denovodb" ~ "green"
        ),
        SHAPE=case_when(
            anno == "LOF"   ~ "1",
            anno == "MIS"   ~ "4"
        )
        ) -> anno_pos_id_study_SAHPE_COLOR
    anno_pos_id_study_SAHPE_COLOR %>% mutate(SHAPE=as.character(SHAPE)) %>% 
        ggplot(aes(x=pos, y=HEIGHT)) + geom_point(aes(shape=anno, color=study),size=2,alpha=0.5) + 
        xlim(1,POS.AA.MAX) + ylim(-1,0) + g + g__theme + 
        theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank()) + 
        scale_color_manual(values=c("DDD"="red","YCU"="blue","denovodb"="green3")) + 
        scale_shape_manual(values=c("LOF"=17,"MIS"=4)) -> g6
    
    ggdraw() + draw_label(ENST) + g + g__theme -> TITLE
    plot_grid(g3, g2, g1, g4, g6, ncol=1, rel_heights=c(0.3,0.5,0.5,0.5,6))
    ggsave(paste0(OUT,"0612.png"),width=4.5,height=8.5,units="cm",dpi=600)
}

extract_sample2 = function(META__DF){
    META__DF %>% mutate(cdna.aa.sample_mis.lof=map_chr(meta,extract_cdna.aa.sample_mis.lof)) %>% select(cdna.aa.sample_mis.lof) %>% unnest %>% separate(cdna.aa.sample_mis.lof,c("anno","pos","id"),";") %>% mutate(pos=as.numeric(pos)) %>% return
}

extract_cdna.aa.sample_mis.lof = function(META__CHAR){ 
    str_split(META__CHAR,";") %>% .[[1]] -> element__vec
    if( length(element__vec)==13 ){ 
        element__vec[13] -> STUDY.SAMPLE
    } else if ( length(element__vec)==14 ){ 
        paste0(element__vec[13],"-",element__vec[14]) -> STUDY.SAMPLE
    } else if ( length(element__vec)==30 ){ 
        paste0("YCU-",element__vec[13]) -> STUDY.SAMPLE
    } else { 
        return(NA)
    }
    element__vec[8]  -> ANNO
    element__vec[9]  -> CDNA
    element__vec[10]  -> AA
    if(AA==""){CDNA -> CDNA.AA}else{AA -> CDNA.AA}
    if(ANNO %in% LOF){"LOF" -> ANNO}else{"MIS" -> ANNO}
    parse_number( str_split(CDNA,"\\.")[[1]][2] ) %/% 3 + 1 -> POS.AA
    paste0(ANNO,";",POS.AA,";",CDNA.AA,"-",STUDY.SAMPLE) %>% return
}

annotate_pfam = function(DF){
    pos.aa__max=c()
    pos.aa__min=c(1)
    pfams=c(DF$pfam[1])
    for(i in 2:nrow(DF)){
        if(DF$pfam[i]!=DF$pfam[i-1]){
            pos.aa__max=c(pos.aa__max, DF$pos.aa[i-1])
            pos.aa__min=c(pos.aa__min, DF$pos.aa[i])
            pfams=c(pfams,DF$pfam[i])
        }
    }
    pos.aa__max=c(pos.aa__max,DF$pos.aa %>% tail(1))
    tibble(pos.aa.start=pos.aa__min,pos.aa.end=pos.aa__max,pfam=pfams) %>% filter(pfam!="nodomain") %>% mutate(y__min=0,y__max=1) %>% return
}

add_pos.aa = function(DF){
    if( DF$strand[[1]] == "+" ){ DF %>% arrange(pos) %>% mutate(n__row=row_number()) -> DF; DF$pos[[1]] -> pos__start; mutate(DF, pos.aa=(n__row-1)%/%3+1 ) %>% return 
    }else if( DF$strand[[1]] == "-" ){ DF %>% arrange(desc(pos)) %>% mutate(n__row=row_number()) -> DF; DF$pos[[1]] -> pos__start; mutate(DF, pos.aa=(abs(n__row-1))%/%3+1 ) %>% return }
}

add_n.exon = function(DF){
    DF %>% arrange(pos.aa) %>% distinct(chr,start,end) %>% mutate(n.exon=row_number()) -> chr_start_end_n.exon
    DF %>% left_join(chr_start_end_n.exon) %>% return
}

MIS=c(
    "missense_variant"
    ,"missense_variant&splice_region_variant"
    ,"protein_protein_contact"
    ,"rare_amino_acid_variant"
)

LOF=c(
    "splice_acceptor_variant&intron_variant"
    ,"splice_acceptor_variant&splice_donor_variant&intron_variant"
    ,"splice_acceptor_variant&splice_region_variant&intron_variant"
    ,"splice_donor_variant&intron_variant"
    ,"splice_donor_variant&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
    ,"exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&3_prime_UTR_variant&intron_variant"
    ,"splice_acceptor_variant&5_prime_UTR_variant&intron_variant"
    ,"splice_acceptor_variant&splice_donor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant"
    ,"splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant"
    ,"splice_acceptor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant"
    ,"splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant"
    ,"stop_gained"
    ,"stop_gained&splice_region_variant"
    ,"start_lost"
    ,"start_lost&splice_region_variant"
    ,"5_prime_UTR_premature_start_codon_gain_variant"
    ,"5_prime_UTR_truncation&exon_loss_variant"
    ,"stop_lost"
    ,"stop_lost&splice_region_variant"
    ,"stop_lost&splice_donor_variant&inframe_deletion&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant"
    ,"start_lost&splice_acceptor_variant&inframe_deletion&splice_region_variant&intron_variant"
    ,"frameshift_variant"
    ,"frameshift_variant&splice_region_variant"
    ,"frameshift_variant&stop_gained"
    ,"frameshift_variant&stop_gained&splice_region_variant"
    ,"frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
    ,"frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant"
    ,"frameshift_variant&splice_donor_variant&intron_variant"
    ,"frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant"
    ,"frameshift_variant&start_lost"
    ,"frameshift_variant&start_lost&splice_region_variant"
    ,"frameshift_variant&stop_lost"
    ,"frameshift_variant&stop_lost&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant"
    ,"frameshift_variant&stop_lost&splice_region_variant"
    ,"splice_donor_variant&inframe_deletion&splice_region_variant&intron_variant"
    ,"splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&splice_donor_variant&inframe_deletion&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
    ,"splice_acceptor_variant&inframe_deletion&splice_region_variant&intron_variant"
    ,"stop_gained&disruptive_inframe_deletion"
    ,"stop_gained&disruptive_inframe_insertion"
    ,"stop_gained&disruptive_inframe_insertion&splice_region_variant"
    ,"stop_gained&inframe_insertion"
    ,"stop_gained&inframe_insertion&splice_region_variant"
    ,"stop_gained&splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
    ,"stop_gained&splice_acceptor_variant&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
    ,"stop_gained&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
)
    
