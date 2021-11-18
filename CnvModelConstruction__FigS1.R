library(tidyverse)
library(cowplot)
gff3_path           ="gencode.v19.annotation.gff3.gz"
protein.enst_path   ="enst.coding.canonical.gencode.noalthaplo.noenstr.txt"
protein.enst.noxy_path="enst.coding.canonical.gencode.noalthaplo.noenstr.noXY.txt"
nc.enst_path        ="enst_nc0130.txt"
nc.enst.noxy_path   ="enst_nc0130.noXY.txt"
sd_path             ="SD.txt"
loeuf_path          ="gnomad.v2.1.1.lof_metrics.by_transcript.txt"
centromere_path     ="centromere.txt"
protein.enst        = read_tsv(protein.enst_path,col_names="enst")
protein.enst.noxy   = read_tsv(protein.enst.noxy_path,col_names="enst")
nc.enst             = read_tsv(nc.enst_path,col_names="enst")
nc.enst.noxy        = read_tsv(nc.enst.noxy_path,col_names="enst")
loeuf               = read_tsv(loeuf_path) %>% select(transcript,oe_lof,oe_lof_lower,oe_lof_upper) %>% rename(enst=transcript)
centromere          = read_tsv(centromere_path)

gff3                            = prep_gencode_gff3(gff3_path)                
gff3.enst                       = prep_gencode_gff3.enst(gff3)               
gff3_ex                         = prep_gencode_gff3_exon(gff3)                
gff3_int                        = prep_gencode_gff3_intron(gff3_ex)           
sd                              = prep_sd(sd_path)                            
covs_list                       = make_covs(gff3,gff3_ex,gff3_int,centromere)    
gff3.enst_enst.len              = covs_list[1]              
gff3.enst_ex.n_len.med_len.sum  = covs_list[2]              
gff3.enst_int.n_len.med_len.sum = covs_list[3]              
gff3.enst_cent.dist_telo.dist   = covs_list[4]
gff3.enst_sd.n                  = make_sd.cov(sd,gff3,out="enst_sd.n_in_R.txt") 
gff3.enst.tolerant              = get_tolerant.region_enst(gff3.enst,loeuf,protein.enst.noxy,margin=100000) 

source("functions__CnvModelConstruction.R")
DATE="xxx"
type="plof"
cnv.n_path			= paste0("dnrate/count_",type,"_enst_mafall.1000k.0309.txt")
out__cor.o.vs.e			= paste0("model_po.col_aic_coef_pval__cor.o.vs.e.",type,DATE,".txt")
out__oe.ratio_snv.vs.cnv	= paste0("model_coef_pval__oe.ratio_snv.vs.cnv.",type,DATE,".txt")
cnv.n_cov_gff3                  = make_cnv.n_covs(cnv.n_path,gff3.enst,correction=TRUE)
cnv.n_cov_protein.noxy          = inner_join(cnv.n_cov_gff3,protein.enst.noxy,by="enst")
cnv.n_cov_protein.noxy.tol      = inner_join(cnv.n_cov_protein.noxy,gff3.enst.tolerant,by="enst")
cnv.n_cov_protein.noxy.intol    = filter(cnv.n_cov_protein.noxy,!enst %in% cnv.n_cov_protein.noxy.tol$enst)
cnv.n_cov_nc.noxy               = filter(cnv.n_cov_gff3,enst %in% nc.enst.noxy$enst)
cnv.n_cov_all.noxy              = rbind(cnv.n_cov_protein.noxy,cnv.n_cov_nc.noxy)
covs                            = setdiff(colnames(cnv.n_cov_protein.noxy.tol),c("cnv.n","enst"))
plot_grid_covs(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy,cnv.n_cov_all.noxy,YMAX=4,memori=TRUE);    ggsave(paste0(DATE,"covs.magni.memori.",type,".1000k.0605.png"),width=18.3,height=9.15,units="cm",dpi=600) 
plot_grid_covs(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy,cnv.n_cov_all.noxy,YMAX=100,memori=TRUE);  ggsave(paste0(DATE,"covs.memori.",      type,".1000k.0605.png"),width=18.3,height=9.15,units="cm",dpi=600)
plot_grid_covs(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy,cnv.n_cov_all.noxy,YMAX=4,memori=FALSE);   ggsave(paste0(DATE,"covs.magni.",       type,".1000k.0605.png"),width=18.3,height=9.15,units="cm",dpi=600)
plot_grid_covs(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy,cnv.n_cov_all.noxy,YMAX=100,memori=FALSE); ggsave(paste0(DATE,"covs.",             type,".1000k.0605.png"),width=18.3,height=9.15,units="cm",dpi=600)

training.data=cnv.n_cov_protein.noxy.tol
test.data=list(cnv.n_cov_protein.noxy.tol,cnv.n_cov_protein.noxy.intol,cnv.n_cov_nc.noxy)
plot_obs.vs.prd_figs(TRAINING=training.data,TESTS=test.data,COLS=c("lightgreen","lightblue","green"),XMAX=3,YMAX=7) -> figs  
plot_oe_snv.vs.cnv_figs(TRAINING=cnv.n_cov_protein.noxy.tol,TEST=cnv.n_cov_protein.noxy,loeuf=loeuf,XMAX=2,YMAX=5)  -> figs2 
plot_grid(figs.plof$fig[[1]],figs.plof$fig[[2]],figs.plof$fig[[3]],figs.cg$fig[[1]],figs.cg$fig[[2]],figs.cg$fig[[3]],nrow=2)
ggsave(paste0(date,"_oe.cor_oe.comp__",type,".png"),width=18.3,height=8.9,units="cm",dpi=600)
plot_grid(figs$fig[[14]],figs$fig[[77]])
ggsave("o.vs.e.png",units="cm",width=9,height=4.5)
plot_grid(figs2$fig[[14]])
ggsave("snv.vs.cnv.png",units="cm",width=4.5,height=4.5)
figs  %>% dplyr::select(-c(test,glm.object,prd,fig)) %>% write_tsv(out__cor.o.vs.e)
figs2 %>% dplyr::select(c(model,coef,pval))     %>% write_tsv(out__oe.ratio_snv.vs.cnv)

all.enst = rbind(protein.enst,nc.enst)
cnv.n_cov_all.enst = inner_join(cnv.n_cov_gff3,all.enst,by="enst")
sapply(c(1:63),function(x){
    get_model(training=cnv.n_cov_protein.noxy,x) %>% add_prd2(TEST=cnv.n_cov_all.enst) %>% write_tsv(paste0("cnv.n_cov_all.enst_prd_",type,"_model",x,"_",DATE,".txt"))
})

