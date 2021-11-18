# plot xhmm score for pc cnvs
library(tidyverse)
library(ggbeeswarm)
library(cowplot)
XhmmSample_path                     = "XhmmSample.txt"
pc_path                             = "pc.result.txt"
nc_path                             = "sample.healthy.txt"
sampleqc_path                       = "sampleqc.20191202.txt"
cnv.id_type_sample_sq_path          = "cnv.id_type_sample_sq.txt"
xcnv_path                           = "DATA.xcnv"
ped_path                            = "invcf.sex.ped"
sample_trio.qc.by.snv_path 	    = "sampleqc.denovo.20191202.txt"
sq.grid                             = c(90) #c(60,70,80,90)
type.grid                           = c("DEL","DUP","TOTAL")
n.target.grid                       = c(1,2,3,4,5)
count.grid                          = c(3)
output.png                          = "sampleqc.cnv.0519.png"
output.sampleqc.txt                 = "sampleqc.cnv.0519.txt"
output.cnv.txt                      = "invcf_cnv.in.xcnv_n.target_sqs.0519.txt"

XhmmSample                          = read_tsv(XhmmSample_path,col_names="sample")
sample_type_cnv.id__for.pc          = read_delim(pc_path,col_types="ccc",col_names=c("sample","type","cnv.id"),delim=" ",comment="#")
sample__for.nc                      = read_tsv(nc_path,col_names=c("sample")) 
cnv.id_type_sample_sq               = read_tsv(cnv.id_type_sample_sq_path,col_names=c("cnv.id","type","sample","sq")) 
xcnv                                = read_tsv(xcnv_path) %>% select(SAMPLE,CNV,Q_SOME,INTERVAL,NUM_TARG) %>% rename(sample=SAMPLE,type=CNV,sq=Q_SOME,cnv.id=INTERVAL,n.target=NUM_TARG)
sample_trio.qc.by.snv      	    = read_tsv(sample_trio.qc.by.snv_path) %>% select(sample,trio.qc.by.snv)
sample_fa_mo                        = read_tsv(ped_path,col_types="_ccc__",col_names=c("sample","fa","mo")) %>% filter(fa!=0) 
read_tsv(sampleqc_path) %>% 
    right_join(XhmmSample) %>% 
    mutate(qced=ifelse(
        coerced.status=="ok" & 
        contami.status=="ok" & 
        sex.sample.status!="outlier" & 
        titvratio.status=="ok" & 
        nsnp.status=="ok" & 
        nins.status=="ok" & 
        ndel.status=="ok" & 
        insdelratio.status=="ok" & 
        hethomratio.status=="ok","ok","outlier")) %>%
    select(sample,qced) %>% rename(qc.by.snv=qced)                          -> XhmmSample_qc.by.snv
get_count.in.control(cnv.id_type_sample_sq,sample__for.nc,SQ.THR=90)        -> cnv.id_type_count.in.control
annotate_cnv.in.trio_with.sq(sample_fa_mo,xcnv,cnv.id_type_sample_sq)       -> sample_fa_mo_type_cnv.id_sq.sample_sq.fa_sq.mo
plot_score(cnv.id_type_sample_sq,sample_type_cnv.id__for.pc,sample__for.nc,sample_fa_mo) -> g1
expand.grid(sq.grid,type.grid,n.target.grid,count.grid) %>% as_tibble %>% 
    rename(sq.thr=Var1,type=Var2,n.target=Var3,count.thr=Var4) %>%
    mutate(data1=list(xcnv),
        data2               = list(sample_fa_mo_type_cnv.id_sq.sample_sq.fa_sq.mo),
        data3               = list(cnv.id_type_count.in.control),
        sample_n.cnv        = pmap(list(SQ.THR=sq.thr,TYPE=type,N.TARGET=n.target,DF=data1),get_cnv.count.table),
        sample_n.dncnv      = pmap(list(SQ.THR=sq.thr,TYPE=type,N.TARGET=n.target,DF=data2),get_denovo.count.table),
        sample_n.rare.dncnv = pmap(list(SQ.THR=sq.thr,TYPE=type,N.TARGET=n.target,COUNT.THR=count.thr,CNV.ID_TYPE_COUNT.IN.CONTROL=data3,DF=data2),get_rare.denovo.count.table)) -> dt

sq.thr                  = 90
n.target                = 4 
n.cnv.del.thr           = 20
n.cnv.dup.thr           = 20 
n.cnv.total.thr         = 30 
n.dncnv.del.thr         = 10 
n.dncnv.dup.thr         = 10 
n.dncnv.total.thr       = 15 
n.rare.dncnv.del.thr    = 3  
n.rare.dncnv.dup.thr    = 3  
n.rare.dncnv.total.thr  = 4  
merge_sampleqc(XhmmSample_qc.by.snv,sample_trio.qc.by.snv,dt,sq.thr,n.target) -> dt2  
write.table(dt2,output.sampleqc.txt,quote=FALSE,row.names=FALSE,sep="\t")
    ggsave(output.png,height=9,width=18,units="cm")
plot_figs(g1,dt2,sq.thr,n.target,output.png)
left_join(sample_fa_mo_type_cnv.id_sq.sample_sq.fa_sq.mo,cnv.id_type_count.in.control) %>% mutate_at(vars("count.in.control"),funs(ifelse(is.na(.),0,.))) %>% write_tsv(output.cnv.txt)

# functions
plot_figs = function(G1,DT2,SQ.THR,N.TARGET,output.png){
    ggplot(DT2,aes(x=n.cnv.del,          fill=qc.by.snv))        + geom_histogram(binwidth=1) + ggtitle(paste("DEL",SQ.THR,N.TARGET)) + xlim(NA,50) -> g1
    ggplot(DT2,aes(x=n.cnv.dup,          fill=qc.by.snv))        + geom_histogram(binwidth=1) + ggtitle(paste("DUP",SQ.THR,N.TARGET)) + xlim(NA,50) -> g2
    ggplot(DT2,aes(x=n.cnv.total,        fill=qc.by.snv))        + geom_histogram(binwidth=1) + ggtitle(paste("TOT",SQ.THR,N.TARGET)) + xlim(NA,50) -> g3
    ggplot(DT2,aes(x=n.dncnv.del,        fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("DEL",SQ.THR,N.TARGET)) + xlim(NA,30) -> g4
    ggplot(DT2,aes(x=n.dncnv.dup,        fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("DUP",SQ.THR,N.TARGET)) + xlim(NA,30) -> g5
    ggplot(DT2,aes(x=n.dncnv.total,      fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("TOT",SQ.THR,N.TARGET)) + xlim(NA,30) -> g6
    ggplot(DT2,aes(x=n.rare.dncnv.del,   fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("DEL",SQ.THR,N.TARGET)) + xlim(NA,10) -> g7
    ggplot(DT2,aes(x=n.rare.dncnv.dup,   fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("DUP",SQ.THR,N.TARGET)) + xlim(NA,10) -> g8
    ggplot(DT2,aes(x=n.rare.dncnv.total, fill=trio.qc.by.snv))   + geom_histogram(binwidth=1) + ggtitle(paste("TOT",SQ.THR,N.TARGET)) + xlim(NA,10) -> g9
    ggplot(DT2,aes(x=n.rare.dncnv.del  ))   + geom_histogram(binwidth=1) + ggtitle(paste("DEL",SQ.THR,N.TARGET)) + xlim(NA,10) -> g7
    ggplot(DT2,aes(x=n.rare.dncnv.dup  ))   + geom_histogram(binwidth=1) + ggtitle(paste("DUP",SQ.THR,N.TARGET)) + xlim(NA,10) -> g8
    ggplot(DT2,aes(x=n.rare.dncnv.total))   + geom_histogram(binwidth=1) + ggtitle(paste("TOT",SQ.THR,N.TARGET)) + xlim(NA,10) -> g9
    plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,nrow=3) -> g1to9
    plot_grid(G1,g1to9,nrow=2)
    plot_grid(G1+g,g1to9,nrow=2,rel_heights=c(2,1)) %>% return
}

merge_sampleqc = function(XhmmSample_qc.by.snv,sample_trio.qc.by.snv,dt,SQ.THR,N.TARGET){
    left_join(XhmmSample_qc.by.snv,sample_trio.qc.by.snv) %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DEL'    & n.target==",N.TARGET)) %>% select(sample_n.cnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DUP'    & n.target==",N.TARGET)) %>% select(sample_n.cnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='TOTAL'  & n.target==",N.TARGET)) %>% select(sample_n.cnv) %>% unnest,by="sample") %>%
        rename(n.cnv.del=n.cnv.x,n.cnv.dup=n.cnv.y,n.cnv.total=n.cnv) %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DEL'    & n.target==",N.TARGET)) %>% select(sample_n.dncnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DUP'    & n.target==",N.TARGET)) %>% select(sample_n.dncnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='TOTAL'  & n.target==",N.TARGET)) %>% select(sample_n.dncnv) %>% unnest,by="sample") %>%
        rename(n.dncnv.del=n.dncnv.x,n.dncnv.dup=n.dncnv.y,n.dncnv.total=n.dncnv) %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DEL'    & n.target==",N.TARGET)) %>% select(sample_n.rare.dncnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='DUP'    & n.target==",N.TARGET)) %>% select(sample_n.rare.dncnv) %>% unnest,by="sample") %>%
        left_join(filter_(dt,paste0("sq.thr==",SQ.THR," & type=='TOTAL'  & n.target==",N.TARGET)) %>% select(sample_n.rare.dncnv) %>% unnest,by="sample") %>%
        rename(n.rare.dncnv.del=n.rare.dncnv.x,n.rare.dncnv.dup=n.rare.dncnv.y,n.rare.dncnv.total=n.rare.dncnv) %>%
        mutate_at(vars(starts_with("n.cnv")),funs(ifelse(is.na(.),0,.))) %>%
        mutate_at(vars(starts_with("n.")),funs(ifelse(!is.na(trio.qc.by.snv) & is.na(.),0,.))) %>%
        mutate(
            n.cnv.del.status            =if_else(n.cnv.del          >n.cnv.del.thr  ,'outlier','ok','NA'),
            n.cnv.dup.status            =if_else(n.cnv.dup          >n.cnv.dup.thr  ,'outlier','ok','NA'),
            n.cnv.total.status          =if_else(n.cnv.total        >n.cnv.total.thr,'outlier','ok','NA'),
            n.dncnv.del.status          =if_else(n.dncnv.del        >n.dncnv.del.thr  ,'outlier','ok','NA'),
            n.dncnv.dup.status          =if_else(n.dncnv.dup        >n.dncnv.dup.thr  ,'outlier','ok','NA'),
            n.dncnv.total.status        =if_else(n.dncnv.total      >n.dncnv.total.thr,'outlier','ok','NA'),
            n.rare.dncnv.del.status     =if_else(n.rare.dncnv.del   >n.rare.dncnv.del.thr  ,'outlier','ok','NA'),
            n.rare.dncnv.dup.status     =if_else(n.rare.dncnv.dup   >n.rare.dncnv.dup.thr  ,'outlier','ok','NA'),
            n.rare.dncnv.total.status   =if_else(n.rare.dncnv.total >n.rare.dncnv.total.thr,'outlier','ok','NA')) %>%
        return
}

annotate_cnv.in.trio_with.sq = function(sample_fa_mo,xcnv,cnv.id_type_sample_sq){
    left_join(sample_fa_mo,xcnv) %>% 
        left_join(cnv.id_type_sample_sq,by=c("cnv.id"="cnv.id","type"="type","fa"="sample")) %>%
        left_join(cnv.id_type_sample_sq,by=c("cnv.id"="cnv.id","type"="type","mo"="sample")) %>%
        rename(sq.sample=sq.x,sq.fa=sq.y,sq.mo=sq) %>%
        mutate_at(vars(starts_with("sq")),funs(ifelse(is.na(.),0,.))) %>% return
}

get_cnv.denovo.rare = function(SQ.THR,TYPE,N.TARGET,cnv.id_type__for.rare.in.control,DF){
    if(TYPE!="TOTAL"){filter_(DF,paste0("type=='",TYPE,"'")) -> DF}
    filter_(DF,paste0("sq.sample>=",SQ.THR," & sq.fa<",SQ.THR," & sq.mo<",SQ.THR," & n.target>=",N.TARGET)) %>% 
        inner_join(cnv.id_type__for.rare.in.control,by=c("type","cnv.id")) %>% 
        full_join(sample.trio_qced) %>% return
}

get_count.in.control = function(CNV.ID_TYPE_SAMPLE_SQ,SAMPLE__FOR.NC,SQ.THR=90){
    inner_join(CNV.ID_TYPE_SAMPLE_SQ,SAMPLE__FOR.NC) %>% 
        mutate(call=if_else(sq>=SQ.THR,1,0)) %>% 
        group_by(cnv.id,type) %>% 
        summarise(count.in.control=sum(call)) %>% return
}

get_rare.denovo.count.table = function(SQ.THR,TYPE,N.TARGET,COUNT.THR,CNV.ID_TYPE_COUNT.IN.CONTROL,DF){
    if(TYPE!="TOTAL"){filter_(DF,paste0("type=='",TYPE,"'")) -> DF}
    filter_(DF,paste0("sq.sample>=",SQ.THR," & sq.fa<",SQ.THR," & sq.mo<",SQ.THR," & n.target>=",N.TARGET)) %>% 
        inner_join(CNV.ID_TYPE_COUNT.IN.CONTROL,by=c("type","cnv.id")) %>% 
        filter(count.in.control<COUNT.THR) %>% group_by(sample) %>% summarise(n.rare.dncnv=n()) %>% return
}

get_denovo.count.table = function(SQ.THR,TYPE,N.TARGET,DF){
    if(TYPE!="TOTAL"){filter_(DF,paste0("type=='",TYPE,"'")) -> DF}
    filter_(DF,paste0("sq.sample>=",SQ.THR," & sq.fa<",SQ.THR," & sq.mo<",SQ.THR," & n.target>=",N.TARGET)) %>%
        group_by(sample) %>% summarise(n.dncnv=n()) %>% return
}

get_cnv.count.table = function(SQ.THR,TYPE,N.TARGET,DF){
    if(TYPE!="TOTAL"){filter_(DF,paste0("type=='",TYPE,"'")) -> DF}
    filter_(DF,paste0("sq>=",SQ.THR," & n.target>=",N.TARGET)) %>% 
        group_by(sample) %>% summarise(n.cnv=n()) %>% return
}

plot_score = function(cnv.id_type_sample_sq,sample_type_cnv.id__for.pc,SAMPLE__FOR.NC,sample_fa_mo){
    inner_join(cnv.id_type_sample_sq,sample_type_cnv.id__for.pc,by=c("cnv.id","type","sample")) -> sample_type_cnv.id_sq__for.pc
    SAMPLE__FOR.NC$sample -> sample.nc.vector
    paste(sample_type_cnv.id__for.pc$cnv.id,sample_type_cnv.id__for.pc$type,sep=";") -> cnv.id.type.vector
    sample_type_cnv.id__for.pc %>% inner_join(sample_fa_mo) %>% select(c(fa,mo)) %>% gather %>% pull(value) -> sample__fa.mo.of.pc
    expand.grid(sample=sample.nc.vector,cnv.id.type.vector) %>% 
        as_tibble %>% 
        separate(Var2,c("cnv.id","type"),";") %>%
        left_join(cnv.id_type_sample_sq) %>% 
        mutate_at(vars(sq),funs(ifelse(is.na(.),0,.))) %>%
        filter(!sample %in% sample__fa.mo.of.pc) -> sample_type_cnv.id_sq__for.nc
    mutate(sample_type_cnv.id_sq__for.pc,order=cnv.id) %>% 
        separate(cnv.id,c("chr","start","end"),convert=T) %>% 
        mutate(size=end-start) %>% 
        distinct(order,type,.keep_all=TRUE) %>% 
        arrange(type,size) %>% 
        pull(order) -> cnv.order
    cnv.order = c(head(cnv.order,n=44) %>% sort,  tail(cnv.order,n=32) %>% sort)
    ggplot(NULL) +
        geom_quasirandom(data=sample_type_cnv.id_sq__for.nc,aes(x=cnv.id,y=sq,color=as.factor(type)),alpha=0.1) +
        geom_point(data=sample_type_cnv.id_sq__for.pc,aes(x=cnv.id,y=sq,color="black"),shape=4) +
        scale_x_discrete(limits=cnv.order) +
        geom_hline(yintercept=90,linetype="dashed") +
        theme(axis.text.x = element_text(angle=90,hjust=1))
}

count_transmission = function(df,SQ.THR){
    filter(df,sq.sample>SQ.THR) %>% mutate(tr=if_else(sq.pt>SQ.THR,"transmitted","not")) %>% count(tr) %>% spread(tr,n) %>% 
        mutate(sq=SQ.THR) %>% return
}

get_highest.sq.transmitted.doubleton = function(cnv.id_type_sample_sq,sample_fa_mo,sample.trio_qced,xcnv,sq.diff.thr=10){
    cnv.id_type_sample_sq %>% nest(-c(cnv.id,type)) %>% 
        mutate(data.sorted=map(data,~arrange(.,desc(sq))),sq.2nd=map_dbl(data.sorted,~.$sq[2]),sq.3rd=map_dbl(data.sorted,~.$sq[3])) %>% 
        mutate_at(vars(starts_with("sq")),funs(ifelse(is.na(.),0,.))) %>% select(-c(data,data.sorted)) -> cnv.id_type_sq.2nd_sq.3rd
    inner_join(sample_fa_mo,sample.trio_qced) %>% rename(pt=sample) %>% mutate(fid=pt) %>% gather(key=relation,value=sample,pt,fa,mo) %>% 
        inner_join(xcnv) %>% nest(-c(cnv.id,type)) %>% 
        inner_join(cnv.id_type_sq.2nd_sq.3rd,by=c("cnv.id","type")) %>%
        mutate(doubleton=map_chr(data,~if_else(nrow(.)==2,"ok","no")),
            one.family=map_chr(data,~if_else( length(unique(.$fid))==1,"ok","no" )),
            mt.sq.3rd=map2_chr(data,sq.3rd,~if_else( min(.x$sq)>.y,"ok","no" ))) %>%
        filter(doubleton=="ok" & one.family=="ok" & mt.sq.3rd=="ok" & sq.2nd-sq.3rd>=sq.diff.thr) %>% 
        unnest %>% select(type,cnv.id,sample) %>% return
}
