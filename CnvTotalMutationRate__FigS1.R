library(tidyverse)
plof_path="sv.name_plof.enst.1000k.0309.txt"
cg_path="sv.name_cg.enst.1000k.0309.txt"
type_asd_enst_path="asd_enst.txt"
enst_path="enst.coding.canonical.gencode.noalthaplo.noenstr.noXY.txt"
enst_prd.plof__path="../cnv.n_cov_all.enst_prd_plof_model14_20200605.txt"
enst_prd.cg__path="../cnv.n_cov_all.enst_prd_cg_model20_20200312.txt"
enst_tol_path="../enst_tol.reg_tol_in_R.1000k.0309.txt"
out.factor="plof.factor_cg.factor.txt"
myenst = read_tsv(enst_path,col_names="enst")
myenst.tolerant = read_tsv(enst_tol_path) %>% inner_join(myenst) 
read_tsv(enst_prd.plof__path,col_types="c___________d") %>% inner_join(myenst.tolerant) %>% pull(prd) %>% sum       -> total.prd.plof
read_tsv(enst_prd.cg__path,col_types="c___________d") %>% inner_join(myenst.tolerant) %>% pull(prd) %>% sum         -> total.prd.cg
mean.ci %>% filter(SVTYPE=="PLOF") %>% pull(mean_ci) %>% .[[1]] %>% filter(bin.af=="singleton") %>% pull(me)        -> n.gene.per.plof.cnv
mean.ci %>% filter(SVTYPE=="CG")   %>% pull(mean_ci) %>% .[[1]] %>% filter(bin.af=="singleton") %>% pull(me)        -> n.gene.per.cg.cnv
(gnomad %>% filter(enst=="tol" & type=="PLOF" & who=="watterson") %>% pull(mean) ) * n.gene.per.plof.cnv / total.prd.plof -> plof.factor 
(gnomad %>% filter(enst=="tol" & type=="CG" & who=="watterson") %>% pull(mean) ) * n.gene.per.cg.cnv / total.prd.cg       -> cg.factor
cbind(plof.factor,cg.factor) %>% as_tibble %>% write_tsv(out.factor)
n.fam=519
svname_type_name_enst = read_tsv(type_asd_enst_path,col_names=c("svname","type","name","enst")) %>% mutate(type=if_else(type=="DEL","PLOF","CG"))
asd.cnv.prop = svname_type_name_enst %>% mutate(
        who=if_else(grepl("p",name),"pt","sib"),
        myenst=map_chr(.$enst,~extract_myenst(.x,myenst$enst)),
        tol.TF=map_lgl(myenst,~check_ensts_are_tolerant(.x,myenst.tolerant$enst,AllOrAny="any")),
        myenst.n = map_dbl(myenst,~strsplit(.,";") %>% unlist %>% length)) %>%
    nest(-c(type,who)) %>% # -> tmp1
    merge_pt_sib %>%
    mutate(
        all=map_dbl(data, ~ filter(.,myenst.n>0) %>% nrow),
        tol=map_dbl(data, ~ filter(.,tol.TF==TRUE) %>% nrow)) %>%
    select(-data) %>%
    gather(key=enst,value=count,all,tol) %>%
    mutate(type_enst=paste(type,enst,sep="_"),
        n.indiv=ifelse(who %in% c("pt","sib"),n.fam,2*n.fam),
	upper=as.numeric(map2(count,n.indiv,~ci95(.x,.y,"up"))),
        lower=as.numeric(map2(count,n.indiv,~ci95(.x,.y,"lo"))),
        mean=count/n.indiv) %>%
    select(-c(count,n.indiv))


order = c("4thmodel","PLOF_pt","PLOF_sib")
out = "tmp20200809.png"
tmp1 %>% mutate(
        data2=map(data, ~prep_boot(.x,n.fam)), 
        mean =map_dbl(data2,mean), 
        up   =map_dbl(data2,~bs_ci(.x,ci="up")), 
        lo   =map_dbl(data2,~bs_ci(.x,ci="lo")),
        type_who=paste0(type,"_",who)) %>% 
    rbind(tibble(type="",who="",data="",data2="",mean=0.01399325,up=0.01399325,lo=0.01399325,type_who="4thmodel")) %>%
    ggplot(aes(x=type_who, y=mean)) +
    geom_point() +
    geom_errorbar(aes(x=type_who, ymax=up, ymin=lo),width=0.01) +
    scale_x_discrete(limits=order) + 
    labs(y="dnCNV rate") + 
    theme(legend.position="none") + 
    g + ylim(0,0.06) 
    ggsave(out,width=4.5,height=4.5,units="cm",dpi=600)

sv_path="gnomad_v2_sv.sites.bed.gz"
enst_path="enst.coding.canonical.gencode.noalthaplo.noenstr.noXY.txt"
out.watterson="20200607watterson.png"
out.n.overlapped.training="20200607n.overlap.training.png"
enst = read_tsv(enst_path,col_names="enst")
plof = read_tsv(plof_path,col_names=c("NAME","plof_gene"))
cg = read_tsv(cg_path,col_names=c("NAME","cg_gene"))
sv <- read_tsv(sv_path, col_types=cols(.default="c")) %>% 
    rename(chrom='#CHROM') %>%
    mutate(SVTYPE=case_when(
        .$NAME %in% plof$NAME ~ "PLOF",
        .$NAME %in% cg$NAME ~ "CG",
        TRUE ~ .$SVTYPE
    )) %>%
    left_join(plof,by="NAME") %>%
    left_join(cg,by="NAME") %>%
    filter(FILTER=="PASS") %>%
    mutate_at(vars(ends_with("_AF")), funs(as.numeric)) %>%
    mutate_at(vars(ends_with("_AC")), funs(as.numeric)) %>%
    mutate_at(vars(ends_with("_AN")), funs(as.numeric)) %>%
    mutate_at(vars(ends_with("_N_BI_GENOS")), funs(as.numeric)) %>%
    mutate_at(vars("AF"), funs(as.numeric)) %>%
    mutate_at(vars("AC"), funs(as.numeric)) %>%
    mutate_at(vars("AN"), funs(as.numeric)) %>%
    mutate_at(vars(contains("FREQ_")), funs(as.numeric)) 
    write_rds(sv,"sv_in_R2.rds")

sv_myenst_tol.TF = sv %>% 
    mutate(myenst_plof=map_chr(.$plof_gene,~extract_myenst(.x,myenst$enst))) %>%
    mutate(myenst_cg=map_chr(.$cg_gene,~extract_myenst(.x,myenst$enst))) %>%
    mutate(myenst=paste(myenst_plof,myenst_cg,sep="")) %>%
    mutate(myenst.n = map_dbl(myenst,~strsplit(.,";") %>% unlist %>% length)) %>%
    mutate(myenst.tol = map_chr(.$myenst,~extract_myenst(.x,myenst.tolerant$enst))) %>%
    mutate(myenst.tol.n = map_dbl(myenst.tol,~strsplit(.,";") %>% unlist %>% length)) %>%
    mutate(tol.TF=map_lgl(myenst,~check_ensts_are_tolerant(.x,myenst.tolerant$enst,AllOrAny="any"))
)

allenst_tolenst_list = list(
    all=filter(sv_myenst_tol.TF,myenst.n>0),
    tol=filter(sv_myenst_tol.TF,tol.TF==TRUE)) %>%
    map(function(x)getAllMus(x,10000))

gnomad = rbind(
    as.data.frame(allenst_tolenst_list$all) %>% rownames_to_column("condition") %>% mutate(enst="all"),
    as.data.frame(allenst_tolenst_list$tol) %>% rownames_to_column("condition") %>% mutate(enst="tol")) %>%
    gather(key=type,value=prop,PLOF,CG) %>%
    mutate(type_enst=paste(type,enst,sep="_"),who="watterson") %>%
    spread(key=condition,value=prop)

order = c("PLOF_tol_watterson","PLOF_tol_all","PLOF_tol_pt","PLOF_tol_sib")

rbind(asd.cnv.prop,gnomad) %>% mutate(type_enst_who=paste0(type_enst,"_",who)) %>% 
    ggplot(aes(x=type_enst_who,y=mean)) +
    geom_point() +
    geom_errorbar(aes(x=type_enst_who,ymax=upper,ymin=lower),width=0.01) +
    scale_x_discrete(limits=order) + 
    labs(y="dnCNV rate") + 
    theme(legend.position="none") + 
    g + ylim(0,0.011) 
    ggsave(out.watterson,width=4.5,height=4.5,units="cm",dpi=600)

library(boot)
library(simpleboot)
mean.ci = mutate(sv_myenst_tol.TF,bin.af=case_when(
        AC==1 ~ "singleton",
        AC==2 ~ "doubleton",
        AF<0.001 ~ "<10E-3",
        AF>=0.001 ~ ">10E-3")) %>%
    filter(tol.TF==TRUE) %>%
    nest(-SVTYPE) %>%
    mutate(
        mean_ci=map(data, ~ group_by(.,bin.af) %>% summarise(me=mean(myenst.tol.n),up=bs_ci(myenst.tol.n,ci="up"),lo=bs_ci(myenst.tol.n,ci="lo")))) %>%
    mutate(mean_ci=map2(SVTYPE,mean_ci,~mutate(.y,type=.x)))

order = c(">10E-3","<10E-3","doubleton","singleton")
rbind(mean.ci$mean_ci[[1]],mean.ci$mean_ci[[2]]) %>%
    ggplot(aes(x=bin.af,y=me,color=type)) +
    geom_point() +
    geom_errorbar(aes(x=bin.af,ymax=up,ymin=lo),width=0.01,position=position_dodge(width=.4)) +
    scale_x_discrete(limits=order) + 
    labs(x="MAF",y="Affected gene N") + 
    theme(legend.position=c(1,1),legend.justification=c(1,1)) + 
    g + ylim(0.5,2.0) 
    ggsave(out.n.overlapped.training,width=4.45,height=4.5,units="cm",dpi=600)


# functions
prep_boot = function(DF, N.FAM){
    DF %>% pull(myenst.n) -> myenst.n__vec
    c( rep(0, N.FAM - length(myenst.n__vec)), myenst.n__vec ) %>% return
}

merge_pt_sib = function(DF){
    filter(DF,type=="PLOF" & who=="pt" ) %>% pull(data) %>% .[[1]] -> plof_pt
    filter(DF,type=="PLOF" & who=="sib") %>% pull(data) %>% .[[1]] -> plof_sib
    filter(DF,type=="CG"   & who=="pt" ) %>% pull(data) %>% .[[1]] -> cg_pt
    filter(DF,type=="CG"   & who=="sib") %>% pull(data) %>% .[[1]] -> cg_sib
    tibble(type=c("PLOF","CG"), who=c("all","all"), data=list(rbind(plof_pt,plof_sib), rbind(cg_pt,cg_sib))) -> DF__new
    rbind(DF, DF__new) %>% return
}

check_ensts_are_tolerant = function(affected_ensts,enst_tolerant_vec,AllOrAny="any"){
    affected_enst_vec =  unlist(strsplit(affected_ensts,";")) 
    tmp = affected_enst_vec %in% enst_tolerant_vec
    if(length(tmp)==0){
        return(FALSE)
    }else{
        if(AllOrAny=="all"){res = all(affected_enst_vec %in% enst_tolerant_vec); return(res)
        }else if(AllOrAny=="any"){res = any(affected_enst_vec %in% enst_tolerant_vec); return(res)
        }    
    }
}

extract_myenst = function(affected_ensts,myenst){
    affected_enst_vec =  unlist(strsplit(affected_ensts,";"))
    affected_myenst_vec = affected_enst_vec[affected_enst_vec %in% myenst]
    return(paste(affected_myenst_vec,collapse=";"))
}

getMu <- function(sv,pop=NULL,svtype=NULL,Ne=10000){
  chr_n = 2*max(select(sv,paste(pop,"_N_BI_GENOS",sep="")),na.rm=T)
  pop_af = paste(pop,'_AF',sep='')
  eval(parse(text=paste0("K=filter(sv,",pop_af,">0 & SVTYPE==svtype & !chrom %in% c('X','Y')) %>% nrow()"))) 
  harmsum <- sum(sapply(1:(chr_n-1),function(k){1/k}))
  theta.hat <- K/harmsum
  mu <- theta.hat/(4*Ne)
  return(mu)
}

getAllMus <- function(sv,Ne){
  svtypes <- c("PLOF","CG")
  mu_pop <- sapply(c("AFR","AMR","EAS","EUR","OTH"),function(pop){
    sapply(svtypes,function(svtype){
      getMu(sv,pop=pop,svtype=svtype,Ne=Ne)
    })
  })
  mu_svtype <- apply(mu_pop,1,function(vals){
    lower <- as.numeric(t.test(vals)$conf.int[1])
    mean <- as.numeric(t.test(vals)$estimate)
    upper <- as.numeric(t.test(vals)$conf.int[2])
    return(c(lower,mean,upper))
  })
  rownames(mu_svtype) <- c("lower","mean","upper")
  return(mu_svtype)
}
 
ci95 = function(count,total,which){
     res = binom.test(count,total)
     if(which=="up"){return(res$conf.int[2])}
     if(which=="lo"){return(res$conf.int[1])}
}

bs_ci = function(df,ci,replicate=5000){
    b.mean = one.boot(df,mean,replicate)
    res = boot.ci(b.mean)
    if(ci=="lo"){return(res$percent[4])}else{return(res$percent[5])}
}

