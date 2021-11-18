library(tidyverse)
library(cowplot)
library(beeswarm)
VCF.PATH="merge.eachvcf.txt"
PATH.VCF_SAMPLE__PATH="path.vcf_sample.txt"
FAM_STATUS_SAMPLE__PATH="nygc_sfari_id_map.csv"
CNV.PC__PATH="pc.dncnv.hg38.txt"
VCF.HEAD=c("chr","start","end","sample","gt","gq","gspc","pl","gp","gl","ft","cn","cnf","cnl","cnp","cnq")
OUT="cnv.id_count_gt.members_denovo.status.txt"
FAM.P1.QC__PATH="fam.p1.qc.permissive.txt"
FAM.S1.QC__PATH="fam.s1.qc.txt"
gq.grid     = seq(0,90,10)
gspc.grid   = 0
count.grid = c(10000)
read_delim(FAM_STATUS_SAMPLE__PATH,delim=",",col_names=c("fam.status","sample"),skip=1) %>% separate(fam.status,c("fam","status")) -> fam_status_sample
read_tsv(PATH.VCF_SAMPLE__PATH,col_names=c("path","sample")) %>% separate(sample,"sample") -> path_sample
inner_join(fam_status_sample,path_sample) %>% nest(-fam) %>% 
    mutate(count=map_dbl(data,function(x){x %>% pull(path) %>% unique %>% length})) %>% 
    mutate(parents=map_lgl(data,function(x){all(c("fa","mo") %in% pull(x,status)) %>% return})) %>% filter(parents==TRUE) %>% { 
        . ->> fam.fa.mo.exist
        unnest(.,data) %>% filter(status=="p1") %>% select(fam) ->> fam.fa.mo.p1.exist
        unnest(.,data) %>% filter(status=="s1") %>% select(fam) ->> fam.fa.mo.s1.exist
        filter(.,count==1) %>% select(fam) ->> fam__in.one.batch 
    }
read_tsv(FAM.P1.QC__PATH,col_types="c") -> fam.p1.qc
read_tsv(FAM.S1.QC__PATH,col_types="c") -> fam.s1.qc
rbind(fam.p1.qc,fam.s1.qc) %>% distinct -> fam.fa.mo.qc 
read_tsv(CNV.PC__PATH,col_type=cols(fam="c")) -> del.pc 

read_tsv(VCF.PATH,col_names=VCF.HEAD,col_types=cols(chr="c")) %>%
    unite("cnv.id",c(chr,start,end),sep="_") %>%
    inner_join(fam_status_sample,by="sample") -> DT  

expand.grid(gq.grid,gspc.grid,count.grid) %>% as_tibble %>%
    dplyr::rename(gq=Var1,gspc=Var2,count.thr=Var3) %>% head(n=10) %>% tail(n=1) %>%
    mutate(
        data1=list(DT),
        data2=pmap(list(DT=data1,GQ=gq,GSPC=gspc),filter_call),
        data2.5=map(data2,~add_count.qced.fa.mo_after_filt(.x,FAM=fam.fa.mo.qc)), 
        data2.7=map2(data2.5,count.thr,~filter(.x,count.cnv<=.y)),
        data3=map(data2.7,spread_gts),  
        data4=map(data3,add_denovo),  
        data5=map(data4,gather_denovo),
        pc.check=map(data5,validate_pc), 
        sensitivity=map_dbl(pc.check,calc_sensitivity),
        denovo.n.p1=map2(data4,"p1",~count_denovo_per_sample(DT=.x,WHO=.y,FAM=fam.fa.mo.p1.exist)),
        denovo.n.s1=map2(data4,"s1",~count_denovo_per_sample(DT=.x,WHO=.y,FAM=fam.fa.mo.s1.exist)), 
        denovo.n.p1.qc=map2(data4,"p1",~count_denovo_per_sample(DT=.x,WHO=.y,FAM=fam.p1.qc)),
        denovo.n.s1.qc=map2(data4,"s1",~count_denovo_per_sample(DT=.x,WHO=.y,FAM=fam.s1.qc)), 
        denovo.mean.p1=map_dbl(denovo.n.p1.qc,mean_denovo), 
        denovo.mean.s1=map_dbl(denovo.n.s1.qc,mean_denovo),
        rate.p1=map2_dbl(data4,list("p1"),~calc_trami(DT=.x,WHO=.y)),
        rate.s1=map2_dbl(data4,list("s1"),~calc_trami(DT=.x,WHO=.y))
        ) -> DT2

filter(DT2,gq==40 & gspc==0 & count.thr==10000) %>% select(denovo.n.p1) %>% pull %>% .[[1]] %>% hist_denovo -> g1 
filter(DT2,gq==40 & gspc==0 & count.thr==10000) %>% select(denovo.n.s1) %>% pull %>% .[[1]] %>% hist_denovo -> g2 
filter(DT2,gq==40 & gspc==0 & count.thr==10000) %>% plot_cnv_size -> g3
plot_1line(DT2,"sensitivity") -> g4
plot_2lines(DT2,"denovo.mean.p1","denovo.mean.s1") -> g5
plot_2lines(DT2,"rate.p1","rate.s1") -> g6
plot_grid(plot_grid(g1+g,g2+g,g3+g,labels=c("a","","b"),nrow=1),plot_grid(g4+g,g5+g+theme(legend.position="none"),g6+g+theme(legend.position="none"),nrow=1),labels=c("","c"),ncol=1)
ggsave("allbatch.permissive.png",width=12.8,height=8.9,units="cm")

filter(DT2,gq==40 & gspc==0 & count.thr==10000) %>% select(data4) %>% unnest %>% 
    filter((denovo.p1=="denovo" & fam %in% fam.p1.qc$fam) | (denovo.s1=="denovo" & fam %in% fam.s1.qc$fam)) %>% 
    write_tsv("cnv.id_count_gt.members_denovo.status__permissive.thr10000.txt")

# functions
plot_beeswarm = function(DT){ggplot(DT,aes(x=gq)) + geom_line(aes(y=VAR1)) %>% return}
plot_1line = function(DT,VAR1){DT %>% ggplot(aes(x=gq)) + geom_line(aes_(y=as.name(VAR1))) %>% return}
plot_2lines = function(DT,VAR1,VAR2){DT %>% gather(key=sample,value=val,VAR1,VAR2) %>% ggplot(aes(x=gq,y=val,color=sample)) + geom_line() %>% return}

plot_cnv_size = function(DT){
    library(ggbeeswarm)
    DT$pc.check[[1]] %>% unnest(pc) %>% mutate(pc.size=end-start) %>% pull(pc.size) -> PC.SIZE.ALL
    DT$pc.check[[1]] %>% unnest(overlap) %>% distinct(chr,start,end) %>% mutate(pc.size=end-start) %>% pull(pc.size) -> PC.SIZE.CALLED
    rbind(cbind(pc.size=setdiff(PC.SIZE.ALL,PC.SIZE.CALLED),type="not.called"),cbind(pc.size=PC.SIZE.CALLED,type="called")) %>%
        as_tibble %>% mutate(pc.size=as.numeric(pc.size)) %>% ggplot(aes(x=type,y=log10(pc.size))) + geom_quasirandom(size=0.5)
}

add_count.qced.fa.mo_after_filt = function(DT,COUNT,FAM){
    DT %>% filter(fam %in% FAM$fam & status %in% c("fa","mo")) %>% 
        mutate(.,gt=case_when(
            str_detect(.$gt, "^0/0") ~ 0,
            str_detect(.$gt, "^0/1") ~ 1,
            str_detect(.$gt, "^1/1") ~ 2, 
            TRUE ~ 0)) %>% group_by(cnv.id) %>% summarise(count.cnv=sum(gt)) -> CNV.ID_COUNT
    left_join(DT,CNV.ID_COUNT) %>% return
}

calc_sensitivity = function(DT){
    DT %>% unnest(overlap) %>% distinct(fam,status,chr,start,end) %>% nrow -> CALLED.PC
    DT %>% nrow -> TOTAL.PC
    CALLED.PC / TOTAL.PC %>% return
}

spread_gts = function(DT){
    select(DT,c(cnv.id,gt,fam,status,count.cnv)) %>% spread(key=status,value=gt) %>% return
}

filter_call = function(DT,GQ,GSPC){
    DT %>% mutate(gt=if_else(gq>=GQ & gspc>=GSPC,gt,"filtered")) %>% return
}

gather_denovo = function(DT){
    DT %>% filter(denovo.p1=="denovo" | denovo.s1=="denovo") %>%
        mutate(
            denovo.p1=ifelse(denovo.p1=="denovo","p1","not"),
            denovo.s1=ifelse(denovo.s1=="denovo","s1","not")) %>%
        gather(key=TMP,value=status,denovo.p1,denovo.s1) %>% 
        filter(status!="not") %>% return
}

validate_pc = function(DT,PC=del.pc){
    PC %>% nest(-c(fam,status),.key="pc") -> PC.nested
    DT %>% select(c(cnv.id,fam,fa,mo,p1,s1,status)) %>% separate(cnv.id,c("chr","start","end"),convert=TRUE) %>% nest(-c(fam,status),.key="calls") -> calls.nested
    left_join(PC.nested,calls.nested) %>%
        mutate(overlap=map2(pc,calls,~detect_overlaped_call(PC2=.x,CALLS2=.y))) %>% return
}

detect_overlaped_call = function(PC2,CALLS2){
    library(GenomicRanges)
    GRanges(seqnames=PC2$chr,ranges=IRanges(PC2$start,end=PC2$end)) -> pc.gr
    GRanges(seqnames=CALLS2$chr,ranges=IRanges(CALLS2$start,end=CALLS2$end)) -> calls.gr
    findOverlaps(pc.gr,calls.gr) -> hits
    pintersect(pc.gr[queryHits(hits)],calls.gr[subjectHits(hits)]) -> overlaps
    width(overlaps)/width(pc.gr[queryHits(hits)]) -> prop.pc
    width(overlaps)/width(calls.gr[subjectHits(hits)]) -> prop.call
    cbind(PC2[queryHits(hits),],CALLS2[subjectHits(hits),],prop.pc,prop.call) %>% as_tibble(.name_repair=make.unique) %>% return
}

count_denovo_per_sample = function(DT,WHO,FAM){
    filter_(DT,paste0("denovo.",WHO,"=='denovo'")) %>% 
        group_by(fam) %>% summarise(count=n()) %>% 
        right_join(FAM) %>% 
        mutate(count=ifelse(is.na(count),0,count)) %>% return
}

hist_denovo = function(DT){
    DT %>% ggplot(aes(x=count)) + geom_histogram(binwidth=1) + ylim(0,25) %>% return
}

mean_denovo = function(DT){
    DT %>% pull(count) %>% mean %>% return
}

add_denovo = function(DT){
    DT %>% mutate(
        denovo.p1=pmap_chr(list(FA=fa,MO=mo,PT=p1),judge_denovo),
        denovo.s1=pmap_chr(list(FA=fa,MO=mo,PT=s1),judge_denovo)) %>% return
}

judge_denovo = function(FA,MO,PT){
    if(is.na(FA) || is.na(MO) || is.na(PT)){return(NA)
        }else if(FA=="filtered" || MO=="filtered" || PT=="filtered"){return(NA)
        }else if(FA=="0/0" && MO=="0/0" && PT=="0/1"){return("denovo")
        }else{return("not")}   
}

calc_trami = function(DT,WHO){
    DT %>% {
        filter_(.,paste('fa=="0/1" & mo=="0/0" & ',WHO,'=="0/1"')) %>% nrow ->> FA.TRANSM
        filter_(.,paste('fa=="0/1" & mo=="0/0" & ',WHO,'=="0/0"')) %>% nrow ->> FA.UNTRAN
        filter_(.,paste('fa=="0/0" & mo=="0/1" & ',WHO,'=="0/1"')) %>% nrow ->> MO.TRANSM
        filter_(.,paste('fa=="0/0" & mo=="0/1" & ',WHO,'=="0/0"')) %>% nrow ->> MO.UNTRAN
    }
            
    (FA.TRANSM+MO.TRANSM)/(FA.TRANSM+FA.UNTRAN+MO.TRANSM+MO.UNTRAN) -> TRANS.PROP            
    return(TRANS.PROP)
}

