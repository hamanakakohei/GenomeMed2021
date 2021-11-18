prep_gencode_gff3 = function(gff3_path){
    gff3 = read_tsv(gff3_path,comment="#",col_types="c_cdd_c_c",col_names=c("chr","type","start","end","strand","att")) %>%
        separate(att,sep=";",c("att1","att2","att3","att4")) %>%
        separate(att4,sep="=",c("n","enst")) %>%    
        separate(enst,sep="\\.",c("enst")) %>%
        mutate(size=end-start+1) %>%
        filter(chr!="chrM") %>% 
        select(-c(att2,att3,n))
    filter(gff3,type=="transcript") %>% group_by(enst) %>% summarise(count=n()) %>% 
        filter(count==2) %>% select(enst) -> enst.par
    filter(gff3,!enst %in% enst.par$enst) -> gff3_no.par
    write.table(gff3_no.par,"gff3_in_R.txt",quote=F,row.names=F,sep="\t")
    return(gff3_no.par)
}

prep_gencode_gff3.enst = function(gff3){
    gff3.enst = filter(gff3,type=="transcript") %>% select(chr,start,end,enst) %>% distinct(enst,.keep_all=TRUE)
    write.table(gff3.enst,"gff3.enst_in_R.txt",quote=F,row.names=F,sep="\t")
    return(gff3.enst)
}

prep_gencode_gff3_exon = function(gff3){
    # input: prep_genecode_gff3 result
    gff3_ex = filter(gff3,type=="exon") %>%
        separate(att1,sep=":",c("n2","n3","ex_th")) %>% 
        mutate(ex_th=as.numeric(ex_th)) %>%
        arrange(enst,ex_th) %>%
        select(-c(n2,n3))
    write.table(gff3_ex,"gff3_ex_in_R.txt",quote=F,row.names=F,sep="\t")
    return(gff3_ex)
}

prep_gencode_gff3_intron = function(gff3_ex){
    # input: prep_genecode_gff3_exon result
    int = data.frame(matrix(rep(NA,6*nrow(gff3_ex)),nrow=nrow(gff3_ex)))
    colnames(int) = c("chr","start","end","size","int_th","enst")
    for(i in 1:(nrow(gff3_ex)-1)){
        if(gff3_ex$enst[i]==gff3_ex$enst[i+1]){
            if(gff3_ex$strand[i]=="+"){
                int$chr[i]=gff3_ex$chr[i]
                int$start[i]=gff3_ex$end[i]+1
                int$end[i]=gff3_ex$start[i+1]-1
                int$size[i]=gff3_ex$start[i+1] - gff3_ex$end[i] - 1
                int$int_th[i]=gff3_ex$ex_th[i]
                int$enst[i]=gff3_ex$enst[i]
            }else{
                int$chr[i]=gff3_ex$chr[i]
                int$start[i]=gff3_ex$end[i+1]+1
                int$end[i]=gff3_ex$start[i]-1
                int$size[i]=gff3_ex$start[i] - gff3_ex$end[i+1] - 1
                int$int_th[i]=gff3_ex$ex_th[i]
                int$enst[i]=gff3_ex$enst[i]
            }
        }
    }
    write.table(int,"gff3_int_in_R.txt",quote=F,row.names=F,sep="\t")
    return(as_tibble(int))
}

prep_sd = function(sd_path){
    # input: segdup file downloaded from UCSC genome browser 
    sd = read_tsv(sd_path,col_types="_cdd__ccdd________________d___") %>% 
        filter(chromStart!=otherStart) %>%  
        mutate(type=case_when(
            .$chrom == .$otherChrom & .$strand == "+" ~ "intra_f",
            .$chrom == .$otherChrom & .$strand == "-" ~ "intra_r",
            .$chrom != .$otherChrom ~ "inter")) %>%
        mutate(
            size_med=if_else(type != "inter", (chromEnd-chromStart+otherEnd-otherStart)/2,NA_real_),
            dist=if_else(type != "inter", ((otherStart+otherEnd)/2)-((chromStart+chromEnd)/2),NA_real_)) %>%
        filter(type!="inter" & chromStart < otherStart) 
    write.table(sd,"sd_in_R.txt",quote=F,row.names=F,sep="\t")
    return(sd)
}

make_covs = function(gff3,gff3_ex,gff3_int,centromere){
    gff3.enst_enst.len = filter(gff3,type=="transcript") %>% 
        select(c(enst,size)) %>%
        dplyr::rename(enst.len=size)
    gff3.enst_ex.n_len.med_len.sum = gff3_ex %>% 
        group_by(enst) %>% 
        summarise(ex.n=n(),ex.len.med=median(size),ex.len.sum=sum(size)) %>% 
        select(c(enst,ex.n,ex.len.med,ex.len.sum))
    gff3.enst_int.n_len.med_len.sum = filter(gff3_int,!is.na(chr)) %>% 
        group_by(enst) %>% 
        summarise(int.n=n(),int.len.med=median(size),int.len.sum=sum(size)) %>% 
        select(c(enst,int.n,int.len.med,int.len.sum))
    gff3.enst_cent.dist_telo.dist = filter(gff3,type=="transcript") %>%
        left_join(centromere,by="chr") %>%
        mutate(pos=(start+end)/2,
            cent.dist=pmap_dbl(list(tel1,int1,cen1,cen2,int2,tel2,pos),calc_cent_dist),
            telo.dist=pmap_dbl(list(tel1,int1,cen1,cen2,int2,tel2,pos),calc_tel_dist),
            cent.dist=if_else(is.na(cent.dist),min(cent.dist,na.rm=T),cent.dist)) %>%
        select(c(enst,cent.dist,telo.dist))
    write.table(gff3.enst_enst.len,"enst_enst.len_in_R.txt",quote=F,row.names=F,sep="\t")
    write.table(gff3.enst_ex.n_len.med_len.sum,"enst_ex.n_len.med_len.sum_in_R.txt",quote=F,row.names=F,sep="\t")
    write.table(gff3.enst_int.n_len.med_len.sum,"enst_int.n_len.med_len.sum_in_R.txt",quote=F,row.names=F,sep="\t")
    write.table(gff3.enst_cent.dist_telo.dist,"enst_cent.dist_telo.dist_in_R.txt",quote=F,row.names=F,sep="\t")
    return(list(gff3.enst_enst.len,gff3.enst_ex.n_len.med_len.sum,gff3.enst_int.n_len.med_len.sum,gff3.enst_cent.dist_telo.dist))
}

calc_cent_dist = function(tel1,int1,cen1,cen2,int2,tel2,pos){
    if(pos<=int1){return(cen1-int1)
    }else if(pos>int1 & pos<=cen1){return(cen1-pos)
    }else if(pos>=cen2 & pos<int2){return(pos-cen2)
    }else if(pos>=int2){return(int2-cen2)
    }else{return(NA)}
}

calc_tel_dist = function(tel1,int1,cen1,cen2,int2,tel2,pos){
    if(pos<=int1){return(pos-tel1)
    }else if(pos>int1 & pos<=cen1){return(int1-tel1)
    }else if(pos>=cen2 & pos<int2){return(tel2-int2)
    }else if(pos>=int2){return(tel2-pos)
    }else{return((tel2-tel1)/2)}
}

make_sd.cov = function(sd,gff3,out){
    library(GenomicRanges)
    sd = dplyr::filter(sd,type=="intra_f")
    sd_gr = GRanges(seqnames=sd$chrom,ranges=IRanges(sd$chromStart,end=sd$otherEnd),score=sd$fracMatch,size=sd$size_med)
    gff3.tr = filter(gff3,type=="transcript")
    gff3.tr_gr = GRanges(seqnames=gff3.tr$chr,IRanges(gff3.tr$start,end=gff3.tr$end,names=gff3.tr$enst))
    gff3.enst_sd.n = dplyr::select(gff3.tr,enst)
    range = c(0,1e7) 
    for(i in 2:length(range)){
        sd_gr_ra = sd_gr[width(sd_gr)>=range[i-1] & width(sd_gr)<range[i]]  
        gff3.enst_sd.n_ra = countOverlaps(gff3.tr_gr,sd_gr_ra) %>% enframe() 
        colnames(gff3.enst_sd.n_ra) = c("enst",paste("sd.n_",as.character(log10(range[i])),sep="")) 
        gff3.enst_sd.n = left_join(gff3.enst_sd.n,gff3.enst_sd.n_ra)
    }
    write.table(gff3.enst_sd.n,out,quote=F,row.names=F,sep="\t")
    return(gff3.enst_sd.n)
}

make_transcript_bed_for_cg = function(gff3,output="transcript_for_cg.bed"){
    filter(gff3,type=="transcript") %>%
        select(c(chr,start,end,enst)) %>%
        mutate(start=start-1) %>%
        mutate(chr=str_replace(chr,pattern="chr",replacement="")) %>%
        write.table(output,col.names=F,row.names=F,quote=F,sep="\t")
}

make_exon_bed_for_plof = function(gff3_ex,output="exon_for_plof.bed"){
    select(gff3_ex,c(chr,start,end,enst)) %>%
        mutate(start=start-1) %>%
        mutate(chr=str_replace(chr,pattern="chr",replacement="")) %>%
        write.table(output,col.names=F,row.names=F,quote=F,sep="\t")
}

get_tolerant.region_enst = function(gff3.enst,loeuf,protein.enst.noxy,margin=1000000){
    inner_join(protein.enst.noxy,loeuf) %>% mutate(bin=as.factor(ntile(oe_lof_upper,10))) %>% group_by(bin) %>% 
        summarise(thr=max(oe_lof_upper)) %>% {
            filter(.,bin==1) %>% select(thr) %>% as.numeric ->> decile12
        }
    gff3.enst_margin = mutate(gff3.enst,sta.mar=start-margin,end.mar=end+margin,sta.mar=if_else(sta.mar<1,1,sta.mar)) 
    select(gff3.enst_margin,c(chr,sta.mar,end.mar,enst)) %>% mutate(chr=str_sub(chr,4,-1)) %>% write_tsv("enst_margin.1000k.0309.bed",quote=F,col_names=F)
    gff3.enst_intol = filter(loeuf,oe_lof_upper<decile12) %>% dplyr::select(enst) %>% inner_join(gff3.enst) 
    
    gff3.enst_margin_gr = GRanges(gff3.enst_margin$chr,IRanges(gff3.enst_margin$sta.mar,gff3.enst_margin$end.mar,names=gff3.enst_margin$enst))
    gff3.enst_intol_gr = GRanges(gff3.enst_intol$chr,IRanges(gff3.enst_intol$start,gff3.enst_intol$end))
    gff3.enst_intol.count = countOverlaps(gff3.enst_margin_gr,gff3.enst_intol_gr)
    names(gff3.enst_intol.count)[gff3.enst_intol.count==0] %>% as_tibble %>% dplyr::rename(enst=value) %>%
        inner_join(filter(loeuf,oe_lof_upper>=decile12)) %>% 
        select(enst) -> gff3.enst_tol.reg_tol
    write.table(gff3.enst_tol.reg_tol,"enst_tol.reg_tol_in_R.1000k.0309.txt",row.names=F,quote=F,sep="\t")
    return(gff3.enst_tol.reg_tol)
}

make_cnv.n_covs = function(cnv.n_path="count_plof_enst_mafall.txt",gff3.enst,correction=TRUE){
    cnv.n = read_delim(cnv.n_path,col_names=c("cnv.n","enst"),delim=" ") %>% 
        mutate(cnv.n=as.integer(str_replace_all(cnv.n,pattern=" ",replacement="")))
    
    cnv.n_cov = dplyr::select(gff3.enst,enst) %>% 
        left_join(cnv.n,by="enst") %>%
        left_join(gff3.enst_enst.len,by="enst") %>%
        left_join(gff3.enst_ex.n_len.med_len.sum,by="enst") %>%
        left_join(gff3.enst_int.n_len.med_len.sum,by="enst") %>%
        left_join(gff3.enst_sd.n,by="enst") %>%
        left_join(gff3.enst_cent.dist_telo.dist,by="enst") %>%
        mutate_at(vars(c(cnv.n,int.n,int.len.med,int.len.sum)),funs(ifelse(is.na(.),0,.)))
    
    if(correction==TRUE){
        cnv.n_cov = mutate(cnv.n_cov,
            enst.len=log(enst.len),
            ex.n=log(ex.n),
            ex.len.sum=log(ex.len.sum),
            cent.dist=log(if_else(cent.dist>1e+7,1e+7,cent.dist)),
            telo.dist=log(if_else(telo.dist>2.5e+7,2.5e+7,telo.dist)),
            sd.n_7=log(sd.n_7+1))
    }
    return(cnv.n_cov)
}

plot_covs = function(cnv.n_cov_tol,covs,n.window=20,n.point=1000,output="0715-2_"){
    for(cov in covs){
        bin_median = return_binned_median(cnv.n_cov_tol,cov,"cnv.n",n.window,n.point)
        ggplot(NULL) +
            eval(parse(text=paste0("geom_point(data=cnv.n_cov_tol,aes(x=",cov,",y=cnv.n),alpha=0.1)"))) +
            geom_line(data=bin_median,aes(x=bin,y=median.each.bin),color="red",size=1) + ylim(0,12) + labs(y="CNV N") + g 
            ggsave(paste0(output,cov,".png"),width=1.75,height=1.75)
    }
}

plot_one_cov = function(cnv.n_cov_tol,cov,n.window=20,n.point=1000,all.data=cnv.n_cov_all.noxy,po.col,YMAX){
    all.data %>% dplyr::select(cov) %>% min(na.rm=T) -> x.min
    return_binned_median(cnv.n_cov_tol,cov,"cnv.n",n.window,n.point) -> bin_median
    return_binned_median(all.data,cov,"cnv.n",n.window,n.point) -> bin_median__all.data
    eval(parse(text=paste0('cor(cnv.n_cov_tol$cnv.n, cnv.n_cov_tol$',cov,' ,use="complete.obs") -> res.corcoeff' )))
    eval(parse(text=paste0('cor.test(cnv.n_cov_tol$cnv.n, cnv.n_cov_tol$',cov,' ) -> res.cortest' )))
    sprintf("%.2f",signif(res.corcoeff,d=2)) -> res.r
    sprintf("%.1e",signif(res.cortest$p.value,d=2)) -> res.p
    bin_median__all.data %>% filter(!is.na(median.each.bin)) %>% pull(bin) %>% max -> x.max
    ggplot(NULL) +
        eval(parse(text=paste0("geom_point(data=cnv.n_cov_tol,aes(x=",cov,",y=cnv.n),alpha=0.1,color=po.col)"))) +
        geom_line(data=bin_median,aes(x=bin,y=median.each.bin),color="black",size=1) + xlim(x.min,x.max) + ylim(0,YMAX) + 
        labs(y="CNV N") + ggtitle(paste0("r: ",res.r,"   p: ",res.p)) + g 
}

return_binned_median = function(df,binned_col,return_col,n.window,n.point,n.point.thr=50){
    lo = dplyr::select(df,binned_col) %>% min(na.rm=T)
    up = dplyr::select(df,binned_col) %>% max(na.rm=T)
    wi = (up-lo)/n.window 
    po = seq(lo,up,length=n.point)
    po_lo = po - wi/2
    po_up = po + wi/2
    medians = sapply(1:n.point,function(i){
            eval(parse(text=paste0("n.point.used=dplyr::filter(df,",binned_col,">po_lo[i] & ",binned_col,"<po_up[i]) %>% nrow")))
            if(n.point.used>n.point.thr){
                eval(parse(text=paste0("median_each_bin=dplyr::filter(df,",binned_col,">po_lo[i] & ",binned_col,"<po_up[i]) %>% summarise(mean(",return_col,"))")))
                return(median_each_bin)
            }else{
                return(NA)
            }
        }
    )
    bin_median.each.bin = as_tibble(cbind(bin=po,median.each.bin=unlist(medians)))
    return(bin_median.each.bin)
}

plot_grid_covs = function(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy,cnv.n_cov_all.noxy,YMAX,memori=TRUE){
    tibble(cov=as.list(covs),data=list(cnv.n_cov_protein.noxy.tol))     %>% mutate(fig=map2(data,cov,~plot_one_cov(.x,.y,all.data=cnv.n_cov_all.noxy,YMAX=YMAX,po.col="red"))) -> figs.tol
    tibble(cov=as.list(covs),data=list(cnv.n_cov_protein.noxy.intol))   %>% mutate(fig=map2(data,cov,~plot_one_cov(.x,.y,all.data=cnv.n_cov_all.noxy,YMAX=YMAX,po.col="lightblue"))) -> figs.intol
    tibble(cov=as.list(covs),data=list(cnv.n_cov_nc.noxy))              %>% mutate(fig=map2(data,cov,~plot_one_cov(.x,.y,all.data=cnv.n_cov_all.noxy,YMAX=YMAX,po.col="green"))) -> figs.nc
    g.no.x = theme(axis.title.x=element_blank(),axis.text.x=element_blank())
    g.no.y = theme(axis.title.y=element_blank(),axis.text.y=element_blank())
    if(memori){
        plot_grid(
            figs.tol$fig[[1]]   ,
            figs.tol$fig[[2]]   ,
            figs.tol$fig[[4]]   ,
            figs.tol$fig[[8]]   ,
            figs.tol$fig[[9]]   ,
            figs.tol$fig[[10]]  ,
            figs.intol$fig[[1]] ,
            figs.intol$fig[[2]] ,
            figs.intol$fig[[4]] ,
            figs.intol$fig[[8]] ,
            figs.intol$fig[[9]] ,
            figs.intol$fig[[10]],
            figs.nc$fig[[1]]    ,
            figs.nc$fig[[2]]    ,
            figs.nc$fig[[4]]    ,
            figs.nc$fig[[8]]    ,
            figs.nc$fig[[9]]    ,
            figs.nc$fig[[10]]   ,
            nrow=3) %>% return
    }else{
        plot_grid(
            figs.tol$fig[[1]]   +g.no.x+g.no.y,
            figs.tol$fig[[2]]   +g.no.x+g.no.y,
            figs.tol$fig[[4]]   +g.no.x+g.no.y,
            figs.tol$fig[[8]]   +g.no.x+g.no.y,
            figs.tol$fig[[9]]   +g.no.x+g.no.y,
            figs.tol$fig[[10]]  +g.no.x+g.no.y,
            figs.intol$fig[[1]] +g.no.x+g.no.y,
            figs.intol$fig[[2]] +g.no.x+g.no.y,
            figs.intol$fig[[4]] +g.no.x+g.no.y,
            figs.intol$fig[[8]] +g.no.x+g.no.y,
            figs.intol$fig[[9]] +g.no.x+g.no.y,
            figs.intol$fig[[10]]+g.no.x+g.no.y,
            figs.nc$fig[[1]]    +g.no.x+g.no.y,
            figs.nc$fig[[2]]    +g.no.x+g.no.y,
            figs.nc$fig[[4]]    +g.no.x+g.no.y,
            figs.nc$fig[[8]]    +g.no.x+g.no.y,
            figs.nc$fig[[9]]    +g.no.x+g.no.y,
            figs.nc$fig[[10]]   +g.no.x+g.no.y,
            nrow=3) %>% return
    }
}

select_cov.comb_low.aic = function(covs,training){
    library(MASS)
    lapply(1:length(covs),function(n){
            combn(x=covs,m=n) %>% as.data.frame %>% as.list %>% lapply(function(comb)paste(comb,collapse="+"))  
        }
    ) %>% unlist -> cov_combs
    names(cov_combs) = as.vector(cov_combs)
    sapply(cov_combs,function(each.comb){
        eval(parse(text = paste0("glm.fit = glm.nb(cnv.n ~ ",each.comb,", data=training)")))
        return(glm.fit$aic)
    }) -> cov.comb_aic
    cov.comb_aic[which.min(cov.comb_aic)] %>% names %>% return
}

get_model = function(training=cnv.n_cov_protein.tol,MODEL){
    library(MASS)
    if(MODEL==1){        glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==2){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==3){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==4){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum          + cent.dist + telo.dist, data=training)
    }else if(MODEL==5){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n +            + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==6){  glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==7){  glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==8){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum + sd.n_7                        , data=training)
    }else if(MODEL==9){  glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum          + cent.dist            , data=training)
    }else if(MODEL==10){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n              + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==11){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==12){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==13){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum                      + telo.dist, data=training)
    }else if(MODEL==14){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n              + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==15){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==16){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==17){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n                       + cent.dist + telo.dist, data=training)
    }else if(MODEL==18){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum          + cent.dist + telo.dist, data=training)
    }else if(MODEL==19){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum          + cent.dist + telo.dist, data=training)
    }else if(MODEL==20){ glm.fit = glm.nb(cnv.n ~ enst.len +                   + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==21){ glm.fit = glm.nb(cnv.n ~            ex.n              + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==22){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum + sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==23){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n + ex.len.sum                                 , data=training)
    }else if(MODEL==24){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n              + sd.n_7                        , data=training)
    }else if(MODEL==25){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum + sd.n_7                        , data=training)
    }else if(MODEL==26){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum + sd.n_7                        , data=training)
    }else if(MODEL==27){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n                       + cent.dist            , data=training)
    }else if(MODEL==28){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum          + cent.dist            , data=training)
    }else if(MODEL==29){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum          + cent.dist            , data=training)
    }else if(MODEL==30){ glm.fit = glm.nb(cnv.n ~ enst.len                     + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==31){ glm.fit = glm.nb(cnv.n ~            ex.n              + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==32){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum + sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==33){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n                                   + telo.dist, data=training)
    }else if(MODEL==34){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum                      + telo.dist, data=training)
    }else if(MODEL==35){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum                      + telo.dist, data=training)
    }else if(MODEL==36){ glm.fit = glm.nb(cnv.n ~ enst.len                     + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==37){ glm.fit = glm.nb(cnv.n ~            ex.n              + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==38){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum + sd.n_7             + telo.dist, data=training)
    }else if(MODEL==39){ glm.fit = glm.nb(cnv.n ~ enst.len                              + cent.dist + telo.dist, data=training)
    }else if(MODEL==40){ glm.fit = glm.nb(cnv.n ~            ex.n                       + cent.dist + telo.dist, data=training)
    }else if(MODEL==41){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum          + cent.dist + telo.dist, data=training)
    }else if(MODEL==42){ glm.fit = glm.nb(cnv.n ~                                sd.n_7 + cent.dist + telo.dist, data=training)
    }else if(MODEL==43){ glm.fit = glm.nb(cnv.n ~ enst.len + ex.n                                              , data=training)
    }else if(MODEL==44){ glm.fit = glm.nb(cnv.n ~ enst.len        + ex.len.sum                                 , data=training)
    }else if(MODEL==45){ glm.fit = glm.nb(cnv.n ~            ex.n + ex.len.sum                                 , data=training)
    }else if(MODEL==46){ glm.fit = glm.nb(cnv.n ~ enst.len                     + sd.n_7                        , data=training)
    }else if(MODEL==47){ glm.fit = glm.nb(cnv.n ~            ex.n              + sd.n_7                        , data=training)
    }else if(MODEL==48){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum + sd.n_7                        , data=training)
    }else if(MODEL==49){ glm.fit = glm.nb(cnv.n ~ enst.len                              + cent.dist            , data=training)
    }else if(MODEL==50){ glm.fit = glm.nb(cnv.n ~            ex.n                       + cent.dist            , data=training)
    }else if(MODEL==51){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum          + cent.dist            , data=training)
    }else if(MODEL==52){ glm.fit = glm.nb(cnv.n ~                                sd.n_7 + cent.dist            , data=training)
    }else if(MODEL==53){ glm.fit = glm.nb(cnv.n ~ enst.len                                          + telo.dist, data=training)
    }else if(MODEL==54){ glm.fit = glm.nb(cnv.n ~            ex.n                                   + telo.dist, data=training)
    }else if(MODEL==55){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum                      + telo.dist, data=training)
    }else if(MODEL==56){ glm.fit = glm.nb(cnv.n ~                                sd.n_7             + telo.dist, data=training)
    }else if(MODEL==57){ glm.fit = glm.nb(cnv.n ~                                       + cent.dist + telo.dist, data=training)
    }else if(MODEL==58){ glm.fit = glm.nb(cnv.n ~ enst.len                                                     , data=training)
    }else if(MODEL==59){ glm.fit = glm.nb(cnv.n ~            ex.n                                              , data=training)
    }else if(MODEL==60){ glm.fit = glm.nb(cnv.n ~                   ex.len.sum                                 , data=training)
    }else if(MODEL==61){ glm.fit = glm.nb(cnv.n ~                                sd.n_7                        , data=training)
    }else if(MODEL==62){ glm.fit = glm.nb(cnv.n ~                                         cent.dist            , data=training)
    }else if(MODEL==63){ glm.fit = glm.nb(cnv.n ~                                                     telo.dist, data=training)
    }
    return(glm.fit)
}

add_prd2 = function(GLM=glm.fit,TEST=cnv.n_cov_protein.tol){
    prd = predict.glm(GLM,newdata=TEST,type="response")
    cnv.n_cov_test_prd = cbind(TEST,prd)
    return(as_tibble(cnv.n_cov_test_prd))
}

plot_obs_prd= function(data=cnv.n_cov_test_prd,po.col="black",XMAX,YMAX,OBJ){
    cor.test(data$cnv.n,data$prd) -> cor.res
    bin_binned.cnv.n = return_binned_median(data,"prd","cnv.n",50,2000)
    colnames(bin_binned.cnv.n) = c("bin","binned.cnv.n")
    if(OBJ=="fig"){
        ggplot(data,aes(y=cnv.n,x=prd)) + 
            geom_point(alpha=0.05,color=po.col) +
            geom_line(data=bin_binned.cnv.n,aes(x=bin,y=binned.cnv.n),color="orange",size=0.5) +
            geom_abline(intercept=0,slope=1,size=0.5,linetype="dashed",color="red") + 
            labs(x="Exp N",y="Obs N") + 
            ggtitle(paste0("r=",signif(cor.res$estimate,d=2)," p=",signif(cor.res$p.value,d=2))) +
            xlim(0,XMAX) + ylim(0,YMAX) + g
    }else if(OBJ=="coef"){
        return(cor.res$estimate)
    }else if(OBJ=="pval"){
        return(cor.res$p.value)
    }
}

plot_obs.vs.prd_figs = function(TRAINING,TESTS,COLS,XMAX,YMAX){
    expand.grid(c(1:63),COLS) %>% as_tibble %>% dplyr::rename(model=Var1,po.col=Var2) %>% mutate(test=case_when(
            .$po.col == "lightgreen"       ~ list(TESTS[[1]]),
            .$po.col == "lightblue" ~ list(TESTS[[2]]),
            .$po.col == "green"     ~ list(TESTS[[3]])
        )) %>%
        mutate(glm.object=map(model,~get_model(training=TRAINING,.x))) %>%
        mutate(prd=map2(glm.object,test,~add_prd2(.x,.y))) %>%
        mutate(aic=map_dbl(glm.object,~.x$aic %>% return)) %>%
        mutate(fig=map2(prd,po.col,~plot_obs_prd(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="fig"))) %>% 
        mutate(coef=map2_dbl(prd,po.col,~plot_obs_prd(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="coef"))) %>% 
        mutate(pval=map2_dbl(prd,po.col,~plot_obs_prd(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="pval"))) %>% 
        mutate(model.details=map(glm.object,spread_glm.object)) %>%
        unnest(model.details) %>%
        return
}

spread_glm.object = function(GLM.OBJECT){
    GLM.OBJECT %>% summary %>% .$coefficients -> TMP
    TMP %>% rownames -> COVARIATES
    
    if("(Intercept)" %in% COVARIATES){ intercept.estimate=TMP["(Intercept)","Estimate"]; intercept.p=TMP["(Intercept)","Pr(>|z|)"]}else{ intercept.estimate=NA; intercept.p=NA}
    if("enst.len"    %in% COVARIATES){  enst.len.estimate=TMP[   "enst.len","Estimate"];  enst.len.p=TMP[   "enst.len","Pr(>|z|)"]}else{  enst.len.estimate=NA;  enst.len.p=NA}
    if("ex.n"        %in% COVARIATES){      ex.n.estimate=TMP[       "ex.n","Estimate"];      ex.n.p=TMP[       "ex.n","Pr(>|z|)"]}else{      ex.n.estimate=NA;      ex.n.p=NA}
    if("ex.len.sum"  %in% COVARIATES){ex.len.sum.estimate=TMP[ "ex.len.sum","Estimate"];ex.len.sum.p=TMP[ "ex.len.sum","Pr(>|z|)"]}else{ex.len.sum.estimate=NA;ex.len.sum.p=NA}
    if("sd.n_7"      %in% COVARIATES){    sd.n_7.estimate=TMP[     "sd.n_7","Estimate"];    sd.n_7.p=TMP[     "sd.n_7","Pr(>|z|)"]}else{    sd.n_7.estimate=NA;    sd.n_7.p=NA}
    if("cent.dist"   %in% COVARIATES){ cent.dist.estimate=TMP[  "cent.dist","Estimate"]; cent.dist.p=TMP[  "cent.dist","Pr(>|z|)"]}else{ cent.dist.estimate=NA; cent.dist.p=NA}
    if("telo.dist"   %in% COVARIATES){ telo.dist.estimate=TMP[  "telo.dist","Estimate"]; telo.dist.p=TMP[  "telo.dist","Pr(>|z|)"]}else{ telo.dist.estimate=NA; telo.dist.p=NA}
        
    tibble(
        intercept.estimate  = intercept.estimate,
        enst.len.estimate   =  enst.len.estimate,
        ex.n.estimate       =      ex.n.estimate,
        ex.len.sum.estimate =ex.len.sum.estimate,
        sd.n_7.estimate     =    sd.n_7.estimate,
        cent.dist.estimate  = cent.dist.estimate,
        telo.dist.estimate  = telo.dist.estimate,
        intercept.p         =   intercept.p,
        enst.len.p          =    enst.len.p,
        ex.n.p              =        ex.n.p,
        ex.len.sum.p        =  ex.len.sum.p,
        sd.n_7.p            =      sd.n_7.p,
        cent.dist.p         =   cent.dist.p,
        telo.dist.p         =   telo.dist.p) %>%
    return
}

plot_oe_snv.vs.cnv_figs = function(TRAINING,TEST,loeuf,XMAX,YMAX){
    tibble(model=c(1:63)) %>%
        mutate(glm.object=map(model,~get_model(training=TRAINING,.x))) %>%
        mutate(prd=map(glm.object,~add_prd2(.x,TEST))) %>%
        mutate(fig=map2(prd,list(loeuf),~compare_snv_cnv_oe2(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="fig"))) %>% 
        mutate(coef=map2_dbl(prd,list(loeuf),~compare_snv_cnv_oe2(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="coef"))) %>% 
        mutate(pval=map2_dbl(prd,list(loeuf),~compare_snv_cnv_oe2(.x,.y,XMAX=XMAX,YMAX=YMAX,OBJ="pval"))) %>% 
        return
}

compare_snv_cnv_oe2 = function(cnv.n_cov_enst_prd,loeuf,XMAX,YMAX,OBJ){
    cnv.n_cov_enst_prd_loeuf = left_join(cnv.n_cov_enst_prd,loeuf,by="enst")
    bin_binned.cnv.n = return_binned_median(cnv.n_cov_enst_prd_loeuf,"oe_lof","cnv.n",50,1000)
    bin_binned.prd   = return_binned_median(cnv.n_cov_enst_prd_loeuf,"oe_lof","prd",50,1000)
    colnames(bin_binned.cnv.n) = c("bin","binned.cnv.n")
    colnames(bin_binned.prd) = c("bin","binned.prd")
    bin_cnv.n_prd = inner_join(bin_binned.cnv.n,bin_binned.prd,by="bin")
    cor.test(cnv.n_cov_enst_prd_loeuf$oe_lof,cnv.n_cov_enst_prd_loeuf$cnv.n/cnv.n_cov_enst_prd_loeuf$prd) -> cor.res 
    sprintf("%.2f",signif(cor.res$estimate,d=2)) -> res.r
    sprintf("%.1e",signif(cor.res$p.value,d=2)) -> res.p
    if(OBJ=="fig"){
        ggplot(NULL) +
            geom_point(data=cnv.n_cov_enst_prd_loeuf,aes(x=oe_lof,y=cnv.n/prd),alpha=0.005) + 
            geom_line(data=bin_cnv.n_prd,aes(x=bin,y=binned.cnv.n/binned.prd),color="orange",size=0.5) +
            geom_abline(intercept=0,slope=1,color="red",linetype="dashed",size=0.5) + 
            labs(x="SNV o/e ratio",y="CNV o/e ratio") + 
            ggtitle(paste0("r: ",res.r,"   p: ",res.p)) + 
            g + xlim(0,XMAX) + ylim(0,YMAX)
    }else if(OBJ=="coef"){
        return(cor.res$estimate)
    }else if(OBJ=="pval"){
        return(cor.res$p.value)
    }
}

