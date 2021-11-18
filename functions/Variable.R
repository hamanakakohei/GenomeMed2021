simplify_denovogear = function(dg_path,output){
    dg <- read_delim(dg_path,delim=" ",col_names=as.character(c(1:45),44),col_types=cols(.default="c")) %>%
        select(seq(3,37,2),40,42,44)
    colnames(dg)=c("CHILD_ID","chr","pos","ref","alt","maxlike_null","pp_null","tgt_null(child/mom/dad)","snpcode","code","maxlike_dnm","pp_dnm","tgt_dnm(child/mom/dad)","lookup","flag","DEPTH_ch","DEPTH_da","DEPTH_mo","MQ_ch","MQ_da","MQ_mo")
    dg %>%
        select(c(CHILD_ID,chr,pos,ref,alt,pp_dnm)) %>%
        write.table(output,row.names=F,col.names=F,quote=F,sep="\t")
}

histogram_titvratio = function(titv_path,output.table,output.plot){
    read_tsv(titv_path,comment="#") %>% select(c(Sample,tiTvRatio)) %>% {
        rename(.,sample=Sample,titvratio=tiTvRatio) %>% write.table(output.table,sep="\t",quote=F,col.names=T,row.names=F)
        ggplot(.,aes(x=tiTvRatio)) + geom_histogram(binwidth=0.01)
    }
    ggsave(output.plot)
}

histogram_n.snv_n.ins_n.del_del.ins.ratio_het.homo.ratio = function(count_path,output.table,output.plot){
    col_to_hist = function(df,width,column){
        ggplot(df,aes(x=count)) + geom_histogram(binwidth=width) + labs(title=column)
    }
    
    library(cowplot)
    widths=c(1,1,1,0.001,0.001)
    columns=c("nSNPs","nInsertions","nDeletions","insertionDeletionRatio","hetHomRatio")
    read_tsv(count_path,comment="#") %>% select(c(Sample,nSNPs,nInsertions,nDeletions,insertionDeletionRatio,hetHomRatio)) -> tmp
    tmp %>% rename(sample=Sample) %>% write.table(output.table,sep="\t",quote=F,col.names=T,row.names=F)
    tmp %>% gather(key=metric,value=count,-Sample) %>% nest(-metric) %>% mutate(fg=pmap(list(df=data,width=as.list(widths),column=as.list(columns)),col_to_hist)) -> dt
    plot_grid(dt$fg[[1]],dt$fg[[2]],dt$fg[[3]],dt$fg[[4]],dt$fg[[5]])
    ggsave(output.plot)
}

modify_consequence_sps = function(df,col_name,MPC.THR=FALSE){
    # 1st arg: dataframe
    # 2nd arg: name of col. for variant consequence
    # return dataframe with modified consequence as "class"
    MIS=c("missense_variant"
        ,"missense_variant&splice_region_variant"
        ,"protein_protein_contact"
        ,"rare_amino_acid_variant"
    )
    SPS=c("splice_acceptor_variant&intron_variant"
        ,"splice_acceptor_variant&splice_donor_variant&intron_variant"
        ,"splice_acceptor_variant&splice_region_variant&intron_variant"
        ,"splice_donor_variant&intron_variant"
        ,"splice_donor_variant&splice_region_variant&intron_variant"
        ,"splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
        ,"exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
    )
    NON=c("stop_gained"
        ,"stop_gained&splice_region_variant"
    )
    SYN=c("stop_retained_variant"
        ,"synonymous_variant"
        ,"initiator_codon_variant"
        ,"splice_region_variant&initiator_codon_variant"
        ,"splice_region_variant&stop_retained_variant"
        ,"splice_region_variant&synonymous_variant"
    )
    INF=c("disruptive_inframe_deletion"
        ,"disruptive_inframe_deletion&splice_region_variant"
        ,"disruptive_inframe_insertion"
        ,"inframe_deletion"
        ,"inframe_insertion"
        ,"disruptive_inframe_insertion&splice_region_variant"
        ,"inframe_deletion&splice_region_variant"
        ,"inframe_insertion&splice_region_variant"
        ,"start_lost&disruptive_inframe_deletion"
        ,"start_lost&disruptive_inframe_insertion"
        ,"start_lost&inframe_deletion"
        ,"start_lost&inframe_deletion&splice_region_variant"
        ,"start_lost&inframe_insertion"
        ,"start_lost&inframe_insertion&splice_region_variant"
        ,"stop_lost&disruptive_inframe_deletion"
        ,"stop_lost&inframe_deletion"
        ,"stop_lost&inframe_deletion&splice_region_variant"
        ,"stop_lost&inframe_insertion"
    )
    FS=c(
        "frameshift_variant"
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
    
    eval(parse(text=paste0("df %>% rename(ano=",col_name,") -> df")))
    df  %>% filter( ano %in% SYN)                   %>% mutate(class="SYN")         -> TMP.syn
    df  %>% filter( ano %in% MIS)                   %>% mutate(class="MIS")         -> TMP.mis
    df  %>% filter( ano %in% MIS & mpc >= 2.0)      %>% mutate(class="MPC2.0.MIS")  -> TMP.mpc2.0.mis
    df  %>% filter( ano %in% NON)                   %>% mutate(class="NON")         -> TMP.non
    df  %>% filter( ano %in% SPS)                   %>% mutate(class="SPS")         -> TMP.sps
    df  %>% filter( ano %in% FS )                   %>% mutate(class="FS")          -> TMP.fs
    df  %>% filter( ano %in% INF)                   %>% mutate(class="INF")         -> TMP.inf
    eval(parse(text=paste("rbind(",paste(ls(,pattern="TMP"),collapse=','),") %>% rename(",col_name,"=ano) %>% return") ))
}
