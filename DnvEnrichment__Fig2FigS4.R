library(tidyverse)
library(cowplot)
sample__FpDnv__path		="sample__FalsePositiveDnv.txt" 
sample__FpDncnv__path		="sample__FalsePositiveDncnv.txt" 
category_requirement__path	="ddg2p.edit.txt"
enst_gene__path			="enst_gene.txt"
sample__target_path		="TargetDiseaseSample__20210223.txt"
sample__undiag.targetdis_path	="sample__undiagnosed.targetdisease20200728.txt" #0529???????
gene_enst__ddg2p.ad.hi__path	="gene_enst__ddg2p.ad.hi.txt"
gene_enst__ddg2p.ad__path	="gene_enst__ddg2p.ad.txt"
clinvar__path			="clinvar_20200608.snpeff.norm.txt"
hgmd__path			="hgmd_pro_2019.3_hg19.snpeff.norm.txt"
maf.ycu__path			="id_maf.ycu.ctrl.txt"
variant_infos__gnomad__path	="burden/variant_infos__gnomad.genome.txt"
dnv.ycu__path			="denovo/denovoqc.20200206.txt"
ENST_OF_INTEREST__path		="enst.coding.canonical.gencode.noalthaplo.noenstr.txt"
chr_pos_depth.v4__path		="DepthV41017M.txt.gz"
chr_pos_depth.v5__path		="DepthV51017M.txt.gz"
chr_pos_depth.v6__path		="DepthV61017M.txt.gz"
chr_start_end_enst__path	="enst_start_end.txt"
exhaustive.variant_annotations__path="exhaustive.vcf.rds"
sampleqc__path			="sampleqc.denovo.20200206.txt"
cnv__path			="ycu_plof_small.bed"
protein.nc.enst_prd.plof__path	="cnv.n_cov_all.enst_prd_plof_model14_20200605.txt"
sampleqc.cnv__path		="sampleqc.cnv.1231.txt"
ddd__path			="ddd.31058case.snpeff2.vcf"
denovodb__path			="denovodb.merge.snpeff.vcf.gz"
study_n__denovodb__path		="study_n.of.probands__denovodb.txt"
study_n__ddd__path		="study_n.of.probands__ddd.txt"
gff3.enst.tol__path		="enst_tol.reg_tol_in_R.txt"
ddg2p.gene_enst__path		="gene_enst__ddg2p.ad.confirmed.txt"
mpc__path			="fordist_constraint_official_mpc_values_v2.txt.gz"
plof.factor__path		="plof.factor_cg.factor.txt"
phenotype__denovodb		=c("autism","congenital_heart_disease","developmentalDisorder","epilepsy","intellectualDisability","sporadic_infantile_spasm_syndrome")
study__excluded.denovodb	=c("DDD_2017","ASD1_2","ASD3","Lelieveld2016")
ssc.study			=c("Iossifov","Krumm","ASD1_2","ASD3","Turner2016","Turner_2017","Werling_2018")
sfari.cnv__path			="sfari_plof.p1.bed"
v4.n        			= 130 
v5.n        			= 772 
v6.n        			= 415
n.ddd       			= 31065 * (5981 / 7700)
n.denovodb  			= 8790 * (1492 / 2179)
vq.snv      			= -29.2
td.snv      			= 5.7 
dn.snv      			= 0.5
dg.snv      			= 0.02
df.snv      			= c("TRUE","FALSE")
vq.indel    			= -7.58
td.indel    			= 5.5
df.indel    			= c("TRUE","FALSE")
out         			= "/betelgeuse01/analysis/hamanaka/png/forrest0225.png"
n__sample.sfari.cnv 		= 2375
n__sample.ycu.cnv 		= 1298

read_tsv(sample__FpDnv__path) 																			    	    -> sample__FpDnv
read_tsv(sample__FpDncnv__path) 																			    -> sample__FpDncnv
read_tsv(sample__target_path, col_names="sample",comment="#")                                                                                                                               -> sample__ycu.target
read_tsv(category_requirement__path) %>% rename(gene=gene_symbol) %>% select(gene,DDD_category, allelic_requirement) %>% distinct %>% nest(-gene) %>% mutate(data=map(data,clean_category_allelic.requirement)) %>% unnest -> gene_category_requirement
read_tsv(enst_gene__path)                                                                                                                                                                   -> enst_gene
read_tsv(sample__undiag.targetdis_path, col_names="sample",comment="#")                                                                                                                     -> sample__undiag.targetdis
read_tsv(gene_enst__ddg2p.ad.hi__path)                                                                                                                                                      -> gene_enst__ddg2p.ad.hi
read_tsv(gene_enst__ddg2p.ad__path)                                                                                                                                                         -> gene_enst__ddg2p.ad
read_tsv(clinvar__path,col_types=cols(chr="c")) %>% mutate(chr=paste0("chr",chr)) %>% select(-mut.rate)                                                                                     -> chr_pos_ref_alt_ano_gene_enst__clinvar
read_tsv(hgmd__path,col_types=cols(chr="c"))    %>% mutate(chr=paste0("chr",chr)) %>% select(-mut.rate)                                                                                     -> chr_pos_ref_alt_ano_gene_enst__hgmd
read_tsv(variant_infos__gnomad__path,col_types=cols(chr="c")) %>% filter(filter=="PASS") %>% select(-filter) %>% rename(non.neuro.AF=non_neuro_AF) %>% mutate(chr=paste0("chr",chr))         -> chr_pos_ref_alt_non.neuro.AF__gnomad
read_tsv(maf.ycu__path,col_names=c("id","maf.ycu.cont")) %>% separate(id,c("chr","pos","ref","alt")) %>% mutate(pos=as.numeric(pos)) %>% mutate(chr=paste0("chr",chr))                       -> chr_pos_ref_alt_maf.ycu.cont__ycu
read_tsv(study_n__denovodb__path,col_names=c("study","n.of.probands"))                                                                                                                       -> study_n.of.probands__denovodb
read_tsv(study_n__ddd__path,col_names=c("study","n.of.probands"))                                                                                                                            -> study_n.of.probands__ddd
read_tsv(ddg2p.gene_enst__path) %>% mutate(ddg2p="yes") %>% select(-gene) %>% distinct                                                                                                      -> enst_ddg2p
read_tsv(chr_pos_depth.v4__path,col_types="cdd",col_names=c("chr","pos","depth.v4"))                                                                                                         -> chr_pos_depth.v4
read_tsv(chr_pos_depth.v5__path,col_types="cdd",col_names=c("chr","pos","depth.v5"))                                                                                                         -> chr_pos_depth.v5
read_tsv(chr_pos_depth.v6__path,col_types="cdd",col_names=c("chr","pos","depth.v6"))                                                                                                         -> chr_pos_depth.v6
read_tsv(ENST_OF_INTEREST__path,col_names="enst")                                                                                                                                            -> enst__of.interest
read_tsv(chr_start_end_enst__path,col_names=c("chr","start","end","enst")) %>% distinct(enst,.keep_all=TRUE)                                                                                 -> chr_start_end_enst 
read_tsv(sampleqc__path) %>% filter(trio.qc.by.snv=="ok" & total.filtered.sampleqc=="ok") %>% select(sample)                                                                                 -> qced.sample 
readRDS(exhaustive.variant_annotations__path) %>% separate(enst,c("enst"),":")                                                                                                               -> chr_pos_ref_alt_mut.rate_ano_gene_enst__exhaustive 
parse_ycu_vcf(dnv.ycu__path,vq.snv,td.snv,dn.snv,dg.snv,df.snv,vq.indel,td.indel,df.indel)                                                                                                   -> dnv_annotations__ycu
parse_ddd_vcf(ddd__path)                                                                                                                                                                     -> dnv_annotations__ddd
parse_denovodb_vcf(denovodb__path,phenotype__denovodb,study__excluded.denovodb,ssc.study)                                                                                                    -> dnv_annotations__denovodb
read_tsv(gff3.enst.tol__path)                                                                                                                                                                -> gff3.enst.tol
read_tsv(mpc__path,col_types="cdcc______________d",skip=1,col_names=c("chr","pos","ref","alt","mpc")) %>% mutate(chr=paste0("chr",chr)) %>% distinct(chr,pos,ref,alt,.keep_all=TRUE)         -> chr_pos_ref_alt_mpc
read_tsv(sampleqc.cnv__path) %>% filter(qc.by.snv=="ok" & trio.qc.by.snv=="ok" & n.cnv.total.status=="ok" & n.dncnv.total.status=="ok" & n.rare.dncnv.total.status=="ok") %>% select(sample) -> sample__dncnv.ycu.qced
read_tsv(cnv__path,col_types="cccc____c",col_names=c("chr","start","end","sample","enst")) %>% mutate(meta=paste(chr,start,end,sep=";")) %>% select(sample,enst,meta) %>% distinct           -> sample_enst_meta__dncnv.ycu 
get_diagnosed_sample(dnv_annotations__ycu, dnv_annotations__ddd, dnv_annotations__denovodb, gene_enst__ddg2p.ad.hi, gene_enst__ddg2p.ad, chr_pos_ref_alt_ano_gene_enst__hgmd, chr_pos_ref_alt_ano_gene_enst__clinvar)  -> sample__diagnosed.ycu.ddd.denovodb
read_tsv(sampleqc__path) %>% filter(trio.qc.by.snv=="ok" & total.filtered.sampleqc=="ok") %>% select(sample)                                                                                 -> sample__ycu.qced
sample__ycu.qced %>% filter(!sample %in% sample__diagnosed.ycu.ddd.denovodb) %>% inner_join(sample__undiag.targetdis)                                                                        -> sample__ycu.qced.undiag 
read_tsv(plof.factor__path) %>% pull(plof.factor)                                                                                                                                            -> plof.factor
read_tsv(protein.nc.enst_prd.plof__path,col_types="c___________d") %>% mutate(prd=prd*plof.factor)                                                                                           -> enst_pred__dncnv
dnv_annotations__ycu      %>% mutate(meta2=meta) %>% separate(meta2,"sample",";") %>% inner_join(sample__ycu.qced) %>% inner_join(sample__ycu.target) %>% anti_join(sample__FpDnv) %>% select(-sample) -> dnv_annotations__ycu    
sample_enst_meta__dncnv.ycu %>% inner_join(sample__dncnv.ycu.qced) %>% inner_join(sample__ycu.qced) %>% inner_join(sample__ycu.target) %>% anti_join(sample__FpDncnv)                        -> sample_enst_meta__dncnv.ycu
n__sample.ycu.cnv + n__sample.sfari.cnv                                                                                                                                                      -> n__dncnv.ycu.sfari

calc_dep.adj.cohort.wide.rate(chr_pos_ref_alt_mut.rate_ano_gene_enst__exhaustive,chr_pos_depth.v4,chr_pos_depth.v5,chr_pos_depth.v6,v4.n,v5.n,v6.n) %>% saveRDS(paste0("prep__ycu.v4.",v4.n,".v5.",v5.n,".v6.",v6.n,".rds")) 
readRDS(paste0("prep__ycu.v4.",v4.n,".v5.",v5.n,".v6.",v6.n,".rds"))                                                                                -> chr_pos_ref_alt_mut.rate_ano_gene_enst_depth.v4_depth.v5_depth.v6_cohort.wide.rate__ycu
calc_cohort.wide.rate(chr_pos_ref_alt_mut.rate_ano_gene_enst__exhaustive,N.TOTAL=n.ddd+n.denovodb)                                                  %>% saveRDS(paste0("prep__ddd.denovodb.",round(n.ddd+n.denovodb),".rds")) 
readRDS(paste0("prep__ddd.denovodb.",round(n.ddd+n.denovodb),".rds"))                                                                                -> chr_pos_ref_alt_mut.rate_ano_gene_enst_cohort.wide.rate__ddd.denovodb 
cohort.wide.rate_at.each.variant_to_at.each.enst.class(chr_pos_ref_alt_mut.rate_ano_gene_enst_depth.v4_depth.v5_depth.v6_cohort.wide.rate__ycu)     %>% saveRDS(paste0("rate__ycu.snv.v4.",v4.n,".v5.",v5.n,".v6.",v6.n,".rds")) 
readRDS(paste0("rate__ycu.snv.v4.",v4.n,".v5.",v5.n,".v6.",v6.n,".rds"))                                                                            -> enst_class_cohort.wide.enst.wide.rate__dnv__ycu 
calc_cohort.wide.enst.wide.dncnv.rate(enst_pred__dncnv,n__dncnv.ycu.sfari)                                                                          %>% saveRDS(paste0("rate__ycu.cnv.",n__dncnv.ycu.sfari,"pt.rds")) 
readRDS(paste0("rate__ycu.cnv.",n__dncnv.ycu.sfari,"pt.rds"))                                                                                       -> enst_class_cohort.wide.enst.wide.rate__dncnv__ycu
rbind(enst_class_cohort.wide.enst.wide.rate__dnv__ycu, enst_class_cohort.wide.enst.wide.rate__dncnv__ycu)                                           %>% saveRDS("rate__ycu.20210223.rds")     
readRDS("rate__ycu.20210223.rds")                                                                                                                   -> enst_class_cohort.wide.enst.wide.rate__ycu
rbind(enst_class_cohort.wide.enst.wide.rate__dnv__ycu)                                                                                              %>% saveRDS("rate__ycu.nocnv.20210223.rds")     
readRDS("rate__ycu.nocnv.0728.rds")                                                                                                                 -> enst_class_cohort.wide.enst.wide.rate__ycu
cohort.wide.rate_at.each.variant_to_at.each.enst.class(chr_pos_ref_alt_mut.rate_ano_gene_enst_cohort.wide.rate__ddd.denovodb)                       %>% saveRDS("rate__ddd.denovodb2.20210223.rds") 
readRDS("rate__ddd.denovodb2.20210223.rds")                                                                                                         -> enst_class_cohort.wide.enst.wide.rate__ddd.denovodb 
cohort.wide.count_at.each.variant_to_at.each.enst.class(dnv_annotations__ycu)                                                                       %>% saveRDS("count__ycu.snv.20210223NoFp.rds") 
readRDS("count__ycu.snv.20210223NoFp.rds")                                                                                                          -> enst_class_cohort.wide.enst.wide.count__dnv__ycu 
calc_cohort.wide.enst.wide.dncnv.count(sample_enst_meta__dncnv.ycu,sample_enst_meta__dncnv.sfari)                                           	    %>% saveRDS("count__ycu.cnv.20210223NoFp.rds") 
readRDS("count__ycu.cnv.20210223NoFp.rds")                                                                                                          -> enst_class_cohort.wide.enst.wide.count__dncnv__ycu
rbind(enst_class_cohort.wide.enst.wide.count__dnv__ycu, enst_class_cohort.wide.enst.wide.count__dncnv__ycu)                                         %>% saveRDS("count__ycu.20210223NoFp.rds")    
readRDS("count__ycu.20210223NoFp.rds")                                                                                                              -> enst_class_cohort.wide.enst.wide.count__ycu
rbind(enst_class_cohort.wide.enst.wide.count__dnv__ycu)                                                                                             %>% saveRDS("count__ycu.nocnv.0728.rds")    
readRDS("count__ycu.nocnv.0728.rds")                                                                                                                -> enst_class_cohort.wide.enst.wide.count__ycu
cohort.wide.count_at.each.variant_to_at.each.enst.class(dnv_annotations__ddd)                                                                       %>% saveRDS("count__ddd.20210223.rds")    
readRDS("count__ddd.20210223.rds")                                                                                                                  -> enst_class_cohort.wide.enst.wide.count__ddd 
cohort.wide.count_at.each.variant_to_at.each.enst.class(dnv_annotations__denovodb)                                                                  %>% saveRDS("count__denovodb.20210223.rds")
readRDS("count__denovodb.20210223.rds")                                                                                                             -> enst_class_cohort.wide.enst.wide.count__denovodb

classes.list=list(
c("MPC2.0.MIS"),
c("SPS", "NON", "FS", "CNV"),
c("SPS", "NON", "FS", "CNV", "MPC2.0.MIS")
)

combine_datasets(enst_class_cohort.wide.enst.wide.count__ycu,enst_class_cohort.wide.enst.wide.count__ddd,"count") %>% combine_datasets(enst_class_cohort.wide.enst.wide.count__denovodb,"count") -> enst_class_cohort.wide.enst.wide.count__ycu.ddd.denovodb
combine_datasets(enst_class_cohort.wide.enst.wide.rate__ycu,enst_class_cohort.wide.enst.wide.rate__ddd.denovodb,"rate") -> enst_class_cohort.wide.enst.wide.rate__ycu.ddd.denovodb
test_enrichment(
    YCU.COUNT   = enst_class_cohort.wide.enst.wide.count__ycu.ddd.denovodb,
    YCU.RATE    = enst_class_cohort.wide.enst.wide.rate__ycu.ddd.denovodb,
    CLASSES.LIST= classes.list,
    ENST__OF.INTEREST= enst__of.interest,
    GENOME.WIDE = FALSE ) -> res 
res %>% saveRDS("res20210302withCNV.rds")
fill_NA_at_known_gene(res, unique(c(gene_enst__ddg2p.ad.hi$enst, gene_enst__ddg2p.ad$enst))) -> res2
make_final_table(res, enst_gene, gene_category_requirement) %>% write_tsv("final.table20210302NoFp.txt")
plot_forrest(class_obs_exp_ci.low_ci.up_p.value,classes.of.interest)

# study qc
CLASSES.LIST=list(c("SYN"),c("MIS"),c("NON"))
read_tsv(sampleqc__path) %>% filter(trio.qc.by.snv=="ok" & total.filtered.sampleqc=="ok") %>% select(sample)                                                                  	-> qced.sample 
parse_ycu_vcf(dnv.ycu__path,vq.snv,td.snv,dn.snv,dg.snv,df.snv,vq.indel,td.indel,df.indel,qced.sample) %>% cohort.wide.count_at.each.variant_to_at.each.enst.class(MAF.THR=1) 	-> enst_class_cohort.wide.enst.wide.count__dnv__ycu 
calc_dep.adj.cohort.wide.rate(chr_pos_ref_alt_mut.rate_ano_gene_enst__exhaustive,chr_pos_depth.v4,chr_pos_depth.v5,chr_pos_depth.v6,v4.n,v5.n,v6.n) %>%                       	saveRDS("prep__ycu20200718.rds") 
readRDS("prep__ycu20200718.rds") %>% cohort.wide.rate_at.each.variant_to_at.each.enst.class(MAF.THR=1)                                                                        	-> enst_class_cohort.wide.enst.wide.rate__dnv__ycu 
test_enrichment(enst_class_cohort.wide.enst.wide.count__dnv__ycu, enst_class_cohort.wide.enst.wide.rate__dnv__ycu, CLASSES.LIST, enst__of.interest, GENOME.WIDE=TRUE)         	-> res.ycu
calc_cohort.wide.rate(chr_pos_ref_alt_mut.rate_ano_gene_enst__exhaustive, N.TOTAL=1) %>% cohort.wide.rate_at.each.variant_to_at.each.enst.class(MAF.THR=1)                    	-> exp.originl
dnv_annotations__ddd      %>% cohort.wide.count_at.each.variant_to_at.each.enst.class(MAF.THR=1)                                                                              	-> enst_class_cohort.wide.enst.wide.count__dnv__ddd
dnv_annotations__denovodb %>% cohort.wide.count_at.each.variant_to_at.each.enst.class(MAF.THR=1)                                                                              	-> enst_class_cohort.wide.enst.wide.count__dnv__denovodb
exp.original %>% mutate(cohort.wide.enst.wide.rate=cohort.wide.enst.wide.rate*31065)                                                                                           	-> enst_class_cohort.wide.enst.wide.rate__dnv__ddd
exp.original %>% mutate(cohort.wide.enst.wide.rate=cohort.wide.enst.wide.rate*8790)                                                                                            	-> enst_class_cohort.wide.enst.wide.rate__dnv__denovodb
test_enrichment(enst_class_cohort.wide.enst.wide.count__dnv__ddd,      enst_class_cohort.wide.enst.wide.rate__dnv__ddd,      CLASSES.LIST, enst__of.interest, GENOME.WIDE=TRUE)	-> res.ddd
test_enrichment(enst_class_cohort.wide.enst.wide.count__dnv__denovodb, enst_class_cohort.wide.enst.wide.rate__dnv__denovodb, CLASSES.LIST, enst__of.interest, GENOME.WIDE=TRUE)	-> res.denovodb
plot_grid( res.ycu$fig[[1]] + g, res.ddd$fig[[1]] + g, res.denovodb$fig[[1]] + g, nrow=1)
ggsave(OUT.DDD,width=14,height=4.5,units="cm")


# functions
fill_NA_at_known_gene = function(DF,ENST__VEC){ DF %>% mutate(p.value = ifelse(enst %in% ENST__VEC, NA, p.value)) }

clean_category_allelic.requirement = function(DF){ tibble( DDD_category = DF$DDD_category %>% paste(collapse=";"), allelic_requirement = DF$allelic_requirement %>% paste(collapse=";") ) %>% return }

plot_two_final.table = function(RES.X,RES.Y,OUT,P.THR=FALSE,XY.THR=FALSE){
    full_join(
            RES.X %>% select(enst,p.bh) %>% rename(p.x=p.bh) %>% distinct,
            RES.Y %>% select(enst,p.bh) %>% rename(p.y=p.bh) %>% distinct ) %>%
        mutate_all(funs(ifelse(is.na(.),1,.))) %>%
        inner_join(enst__of.interest) %>%  
        ggplot(aes(-log10(p.x),y=-log10(p.y))) +
        geom_abline(intercept=0,slope=1,linetype="dashed",alpha=0.25) +
        geom_point(size=1,alpha=0.3) + g -> G
        if(XY.THR){ G + xlim(0,XY.THR) + ylim(0,XY.THR) -> G}
	if(P.THR){G + geom_vline(xintercept=P.THR,linetype="dashed",alpha=0.25) + geom_hline(yintercept=P.THR,linetype="dashed",alpha=0.25) -> G} 	
	G
        ggsave(OUT,width=9,height=9,units="cm",dpi=900)
}

plot_two_res = function(RES.X,RES.Y,P__THR=-log10(8.33e-7),OUT="/betelgeuse01/analysis/hamanaka/png/tmp.compare20200803.png"){
    full_join(
            RES.X %>% select(enst,class,p.value) %>% rename(p.x=p.value),
            RES.Y %>% select(enst,class,p.value) %>% rename(p.y=p.value) ) %>%
        mutate_all(funs(ifelse(is.na(.),1,.))) %>%
        inner_join(enst__of.interest) %>% nest(-enst) %>% mutate(data2=map(data, ~ .x %>% arrange(p.y) %>% head(1)) ) %>% unnest(data2) %>% 
        ggplot(aes(-log10(p.x),y=-log10(p.y))) +
        geom_rect(aes(xmin=0,xmax=P__THR,ymin=P__THR,ymax=Inf),fill="pink",alpha=0.03) +
        geom_vline(xintercept=P__THR,linetype="dashed",alpha=0.25) +
        geom_hline(yintercept=P__THR,linetype="dashed",alpha=0.25) +
        geom_abline(intercept=0,slope=1,linetype="dashed",alpha=0.25) +
        geom_point(size=1,alpha=0.3) + g + xlim(0,12) + ylim(0,12)
        ggsave(OUT,width=4.5,height=4.5,units="cm",dpi=900)
}

make_qqman_compare = function(QQMAN,QQMAN.CNV){
    QQMAN.CNV %>% rename(P.cnv=P) -> QQMAN.CNV
    QQMAN.CNV %>% left_join(QQMAN,by=c("SNP","CHR","BP")) %>% mutate(cnv=ifelse(P.cnv<P & P.cnv<1e-6, "TRUE","FALSE")) -> TMP
    TMP %>% filter(cnv=="TRUE") -> QQMAN_gene.cont.cnv
    QQMAN_gene.cont.cnv %>% mutate(SNP=paste0(SNP," (CNV)")) %>% select(SNP,CHR,BP,P.cnv) %>% rename(P=P.cnv) -> QQMAN_gene.cont.cnv.incl.cnv
    QQMAN_gene.cont.cnv %>% select(SNP,CHR,BP,P) -> QQMAN_gene.cont.cnv.excl.cnv
    TMP %>% filter(cnv=="FALSE") %>% select(SNP,CHR,BP,P.cnv) %>% rename(P=P.cnv) -> QQMAN_gene.nocnv
    rbind(QQMAN_gene.cont.cnv.incl.cnv,QQMAN_gene.cont.cnv.excl.cnv,QQMAN_gene.nocnv) %>%
        filter(!SNP %in% c("FOXD4L6 (CNV)","PSMC3 (CNV)","CBWD6 (CNV)")) %>% return
}


func.manhattan = function(DF,GENE.N,HIGHLIGHT,WIDTH=19,HEIGHT=100,GENENAME=TRUE,OUT){
    par(ps=6)
    PVAL = 0.05/GENE.N
    par(omi=c(0,0,0,0))
    png(OUT, width=WIDTH,height=HEIGHT,units="cm",res=600)
    if(GENENAME){
        manhattan(DF,suggestiveline=FALSE,genomewideline=-log10(PVAL),ylim=c(5,25),cex=1,annotatePval=PVAL,highlight=HIGHLIGHT,annotateTop=FALSE)
    }else{
        manhattan(DF,suggestiveline=FALSE,genomewideline=-log10(PVAL),ylim=c(5,25),cex=1,highlight=HIGHLIGHT,annotateTop=FALSE)
    }
    dev.off()
}


adjust_BF_and_BH = function(DF){ DF %>% mutate(p.bf=ifelse(p.value*n.test>1,1,p.value*n.test),  p.bh=p.adjust(p.bf,"BH")) %>% return }

make_qqman_data = function(RES, ENST__OF.INTEREST, ENST_GENE, CHR_START_END_ENST){
    library(reticulate)
    SCIPY = reticulate:::import(module = "scipy.stats")
    RES %>% filter(exp!=0) %>% inner_join(ENST__OF.INTEREST) %>% arrange(p.value) %>% nest(-enst) %>% mutate(data2=map(data, ~ .x %>% arrange(desc(obs),p.value) %>% head(1)) ) %>% unnest(data2) %>% 
        left_join(ENST_GENE) %>% left_join(CHR_START_END_ENST) %>% filter(!chr %in% c("chrY","chrM")) %>%
        mutate(
            chr=if_else(chr=="chrX","chr23",chr),pos=(start+end)/2,chr=as.numeric(str_replace(chr,pattern=c("chr"),replacement=c(""))),
            gene=case_when(
                enst=="ENST00000498285" ~ "VAMP2",
                enst=="ENST00000222270" ~ "KMT2B",
                enst=="ENST00000354995" ~ "CBWD6",
                enst=="ENST00000376358" ~ "WDR45",
                enst=="ENST00000354995" ~ "CBWD6",
                enst=="ENST00000601123" ~ "AC018630.1",
                enst=="ENST00000570054" ~ "VPS4A",
                enst=="ENST00000538264" ~ "AL592284.1",
                enst=="ENST00000398092" ~ "ARF3",
                enst=="ENST00000335616" ~ "CTD-2600O9.1",
                TRUE ~ gene
            )
            ) %>%
        rename(SNP=gene,CHR=chr,BP=pos,P=p.value) %>% select(obs,exp,SNP,CHR,BP,P) %>% 
        mutate(P=map2_dbl(obs,exp,~SCIPY$poisson$sf(.x-1,.y))) %>% select(-c(obs,exp)) %>% return
}

add_some_gene_name = function(DF){
    DF %>% mutate(
            gene=case_when(
                enst=="ENST00000498285" ~ "VAMP2",
                enst=="ENST00000222270" ~ "KMT2B",
                enst=="ENST00000354995" ~ "CBWD6",
                enst=="ENST00000376358" ~ "WDR45",
                enst=="ENST00000354995" ~ "CBWD6",
                enst=="ENST00000601123" ~ "AC018630.1",
                enst=="ENST00000570054" ~ "VPS4A",
                enst=="ENST00000538264" ~ "AL592284.1",
                enst=="ENST00000398092" ~ "ARF3",
                enst=="ENST00000335616" ~ "CTD-2600O9.1",
                TRUE ~ gene
            )
        ) %>% return
}

make_final_table = function(RES, ENST_GENE, GENE_CATEGORY_REQUIREMENT){
    library(reticulate)
    SCIPY = reticulate:::import(module = "scipy.stats")
    RES %>% filter(exp!=0) %>% mutate(p.value=map2_dbl(obs,exp,~SCIPY$poisson$sf(.x-1,.y))) %>% #filter(obs>0) 
	arrange(p.value) %>% nest(-enst) %>% mutate(n.test=map_dbl(data,nrow), n.test=ifelse(n.test==2,1,n.test), data2=map(data, ~ .x %>% arrange(desc(obs),p.value) %>% head(1)) ) %>% unnest(data2) %>% 
        select(-c(data,poisson.res,ddg2p,ci.low,ci.up)) %>% adjust_BF_and_BH %>% unnest(meta) %>% filter(meta!="NULL") %>% 
        unnest(meta) %>% left_join(ENST_GENE) %>% left_join(GENE_CATEGORY_REQUIREMENT) %>% add_some_gene_name %>% filter(exp!=0) %>% return
}

get_diagnosed_sample = function(YCU,DDD,DENOVODB,DDG2P__HI,DDG2P__AD,HGMD,CLINVAR){
    rbind( 
        select(YCU,      c(chr,pos,ref,alt,ano,enst,meta)) %>% mutate(meta=map_chr(meta,~str_split(.x,";") %>% .[[1]] %>% .[1]  )) , 
        select(DDD,      c(chr,pos,ref,alt,ano,enst,meta))  , 
        select(DENOVODB, c(chr,pos,ref,alt,ano,enst,meta))    ) %>% 
        mutate(mpc=0) %>% 
        modify_consequence_sps("ano",START.TERMINAL=FALSE,UTR=FALSE) %>% {
            filter(., enst %in% DDG2P__HI$enst & class %in% c("NON","SPS","FS")) %>% pull(meta) ->> sample__hi
            filter(., enst %in% DDG2P__AD$enst) %>% inner_join(HGMD, by=c("chr","pos","ref","alt")) %>% pull(meta) ->> sample__ad.hgmd
            filter(., enst %in% DDG2P__AD$enst) %>% inner_join(CLINVAR, by=c("chr","pos","ref","alt")) %>% pull(meta) ->> sample__ad.clinvar
        }
    c(sample__hi, sample__ad.hgmd, sample__ad.clinvar) %>% unique %>% return
}

calc_cohort.wide.enst.wide.dncnv.rate = function(ENST_PRED__DNCNV,N.CNVQC){
    ENST_PRED__DNCNV %>% mutate(cohort.wide.enst.wide.dncnv.rate=prd * N.CNVQC, class="CNV") %>% rename(cohort.wide.enst.wide.rate=cohort.wide.enst.wide.dncnv.rate) %>% select(-prd) %>% return
}

calc_cohort.wide.enst.wide.dncnv.count = function(SAMPLE_ENST__DNCNV,SAMPLE_ENST_META__DNCNV.SFARI=NULL){
    SAMPLE_ENST__DNCNV %>% rbind(SAMPLE_ENST_META__DNCNV.SFARI) %>% group_by(enst) %>% summarise(cohort.wide.enst.wide.dncnv.count=n()) %>% mutate(class="CNV") %>% rename(cohort.wide.enst.wide.count=cohort.wide.enst.wide.dncnv.count) -> enst_cohort.wide.enst.wide.count_class
    SAMPLE_ENST__DNCNV %>% rbind(SAMPLE_ENST_META__DNCNV.SFARI) %>% mutate(meta=paste(sample,meta,sep=";")) %>% select(-sample) %>% nest(.key=meta,-enst) -> enst_meta
    inner_join(enst_cohort.wide.enst.wide.count_class,enst_meta) %>% return
}

cohort.wide.count_at.each.variant_to_at.each.enst.class__simplified = function(DNV_ANNOTATIONS,ANO.COL="ano",STA.STO.TF=FALSE,UTR.TF=FALSE){
    DNV_ANNOTATIONS %>% 
        calc_cohort.wide.enst.wide.dnv.count %>%
        rename(cohort.wide.enst.wide.count=cohort.wide.enst.wide.dnv.count) %>%
        return
}

cohort.wide.count_at.each.variant_to_at.each.enst.class = function(DNV_ANNOTATIONS,MPC=chr_pos_ref_alt_mpc,ANO.COL="ano",STA.STO.TF=FALSE,UTR.TF=FALSE,MPC.THR1=2,MAF.THR=0.00001){
    DNV_ANNOTATIONS %>% 
        left_join(chr_pos_ref_alt_maf.ycu.cont__ycu) %>%
        left_join(chr_pos_ref_alt_non.neuro.AF__gnomad) %>%
        mutate_at(vars(c(maf.ycu.cont,non.neuro.AF)),funs(ifelse(is.na(.),0,.))) %>% 
        filter(maf.ycu.cont<MAF.THR & non.neuro.AF<MAF.THR) %>%
        left_join(MPC) %>% 
        mutate_at(vars(c(maf.ycu.cont,non.neuro.AF,mpc)),funs(ifelse(is.na(.),0,.))) %>% 
        modify_consequence_sps(ANO.COL,START.TERMINAL=STA.STO.TF,UTR=UTR.TF,MPC.THR=MPC.THR1) %>%
        mutate(meta=paste(chr,pos,ref,alt,mpc,ano,variant,enst,gene,meta,sep=";")) %>%
        select(c(enst,class,meta)) %>%
        calc_cohort.wide.enst.wide.dnv.count %>%
        rename(cohort.wide.enst.wide.count=cohort.wide.enst.wide.dnv.count) %>%
        return
}

cohort.wide.rate_at.each.variant_to_at.each.enst.class = function(COHORT.WIDE.RATE_AT.EACH.VARIANT,ANO.COL="ano",STA.STO.TF=FALSE,UTR.TF=FALSE,MPC.THR1=2,MAF.THR=0.00001){
    COHORT.WIDE.RATE_AT.EACH.VARIANT %>%
        left_join(chr_pos_ref_alt_maf.ycu.cont__ycu) %>%
        left_join(chr_pos_ref_alt_non.neuro.AF__gnomad) %>%
        mutate_at(vars(c(maf.ycu.cont,non.neuro.AF)),funs(ifelse(is.na(.),0,.))) %>% 
        filter(maf.ycu.cont<MAF.THR & non.neuro.AF<MAF.THR) %>%
        left_join(chr_pos_ref_alt_mpc) %>% 
        mutate_at(vars(mpc),funs(ifelse(is.na(.),0,.))) %>% 
        modify_consequence_sps(ANO.COL,START.TERMINAL=STA.STO.TF,UTR=UTR.TF,MPC.THR=MPC.THR1) %>%
        add.indel.rate_and_summarise.by.enst.class %>%
        return
}

parse_ycu_vcf = function(DNVYCU__PATH,VQ.SNV,TD.SNV,DN.SNV,DG.SNV,DF.SNV,VQ.INDEL,TD.INDEL,DF.INDEL){
    read_tsv(DNVYCU__PATH,col_types=cols(chr="c")) %>% rename(enst=snpeff.enst,ano=snpeff.csq,gene=snpeff.gene,variant=snpeff.variant) %>% 
        filter_dnv(vq.snv=VQ.SNV,td.snv=TD.SNV,dn.snv=DN.SNV,dg.snv=DG.SNV,df.snv=DF.SNV,vq.indel=VQ.INDEL,td.indel=TD.INDEL,df.indel=DF.INDEL) %>%
        mutate(chr=paste0("chr",chr)) %>%
        mutate(meta=paste(sample,chr,pos,ref,alt,ano,gene,variant,vqslod,triodenovo,dnmfilter,denovogear,denovofilter,SIFT_score,Polyphen2_HVAR_score,CADD_phred,snp20160620_tommo_exome,sep=";")) %>%
        select(c(chr,pos,ref,alt,ano,gene,variant,enst,meta)) %>% return
}

parse_denovodb_vcf = function(DENOVODB__PATH,PHENOTYPE,EXCLUDED.STUDY,SSC.STUDY){
    read_tsv(DENOVODB__PATH,comment="##",col_types="cd_cc__c_ccddccc") %>% rename(chr=`#CHROM`,pos=POS,ref=REF,alt=ALT,study=StudyName,sample=SampleID) %>% 
        separate(INFO,c("n1","ano","n2","gene","n3","n4","enst","n5","n6","cdna","aa"),sep="\\|") %>% select(-c(n1,n2,n3,n4,n5,n6)) %>% mutate(variant=paste(cdna,aa,sep=";")) %>%
        filter(!study %in% EXCLUDED.STUDY & PrimaryPhenotype %in% PHENOTYPE) %>%
        select(c(chr,pos,ref,alt,ano,gene,enst,variant,sample,study)) -> TMP__each.study

    TMP__each.study %>% filter(study %in% SSC.STUDY) %>% 
        distinct(chr,pos,ref,alt,ano,gene,enst,variant,sample) %>% 
        mutate(study="ssc.all") -> TMP__ssc.all
    
    TMP__each.study %>% filter(study %in% c("Yuen2017","Yuen2016")) %>% 
        distinct(chr,pos,ref,alt,ano,gene,enst,variant,sample) %>% 
        mutate(study="yuen.all") -> TMP__yuen.all
    
    TMP__each.study %>% filter(!study %in% SSC.STUDY & !study %in% c("Yuen2017","Yuen2016")) -> TMP__each.study.except.ssc
    
    rbind(TMP__ssc.all, TMP__yuen.all, TMP__each.study.except.ssc) %>%
        mutate(meta=paste(study,sample,sep=";")) %>%
        mutate(chr=paste0("chr",chr)) %>% 
        select(-c(study,sample)) %>%
        return
}

parse_ddd_vcf = function(DDD__PATH){
    read_tsv(DDD__PATH,comment="#",col_types="cdccc__c",col_names=c("chr","pos","study","ref","alt","snpeff")) %>% 
        separate(snpeff,into=c("n1","ano","n2","gene","n3","n4","enst","n5","n6","cdna","aa"),sep="\\|") %>% select(-c(n1,n2,n3,n4,n5,n6)) %>% 
        mutate(variant=paste(cdna,aa,sep=";")) %>%
        select(c(chr,pos,ref,alt,ano,gene,enst,variant,study)) %>% 
        mutate(chr=paste0("chr",chr)) %>%
        rename(meta=study) %>%
        return
}

calc_cohort.wide.rate = function(CHR_POS_REF_ALT_MUT.RATE_ANO_GENE_ENST,N.TOTAL){ CHR_POS_REF_ALT_MUT.RATE_ANO_GENE_ENST %>% mutate(cohort.wide.rate=N.TOTAL*2*mut.rate) %>% return }

calc_dep.adj.cohort.wide.rate = function(CHR_POS_REF_ALT_MUT.RATE_ANO_GENE_ENST,CHR_POS_V4.DEP,CHR_POS_V5.DEP,CHR_POS_V6.DEP,V4.N,V5.N,V6.N){
    CHR_POS_REF_ALT_MUT.RATE_ANO_GENE_ENST %>%
        left_join(CHR_POS_V4.DEP,by=c("chr","pos")) %>%
        left_join(CHR_POS_V5.DEP,by=c("chr","pos")) %>%
        left_join(CHR_POS_V6.DEP,by=c("chr","pos")) %>%
        mutate_at(vars(starts_with("depth.v")),funs(ifelse(is.na(.),0,.))) %>%
        mutate(
            depth.v4=if_else(depth.v4>40,40,depth.v4), #exhaustive ga nan bp aruka, depth ga nan bp aruka, no taiou kankei
            depth.v5=if_else(depth.v5>40,40,depth.v5),
            depth.v6=if_else(depth.v6>40,40,depth.v6)) %>% 
        mutate(dep.adj.cohort.wide.rate=V4.N*2*mut.rate*depth.v4/40+V5.N*2*mut.rate*depth.v5/40+V6.N*2*mut.rate*depth.v6/40) %>%
        rename(cohort.wide.rate=dep.adj.cohort.wide.rate) %>%
        return()
}

add.indel.rate_and_summarise.by.enst.class = function(COHORT.WIDE.RATE_AT.EACH.VARIANT,fs.non.ratio=1.1,fs.inf.ratio=9){
    COHORT.WIDE.RATE_AT.EACH.VARIANT %>%
        group_by(enst,class) %>%
        summarise(cohort.wide.enst.wide.rate=sum(cohort.wide.rate)) %>%
        ungroup %>%
        spread(key=class,value=cohort.wide.enst.wide.rate) %>%
        mutate(FS=fs.non.ratio*NON,INF=FS/fs.inf.ratio) %>% 
        gather(key=class, value=cohort.wide.enst.wide.rate,-enst) %>%
        mutate_at(vars(cohort.wide.enst.wide.rate),funs(if_else(is.na(.),0,.))) %>% 
        return
}

calc_cohort.wide.enst.wide.dnv.count = function(dnvs){
    dnvs %>% nest(.key=meta,-c(enst,class)) %>% mutate(cohort.wide.enst.wide.dnv.count=map_dbl(meta,nrow)) %>% return
}

filter_dnv = function(dnv_annotations2__ycu,vq.snv,td.snv,dn.snv,dg.snv,df.snv,vq.indel,td.indel,df.indel){
    dnv_annotations2__ycu %>%
        left_join(enst_chr_start_end) %>%
        mutate(filter=case_when(
            vqslod>vq.snv & triodenovo>td.snv & dnmfilter>dn.snv & denovogear>dg.snv & denovofilter %in% df.snv & type=="snv" ~ "PASS",
            vqslod>vq.indel & triodenovo>td.indel & denovofilter %in% df.indel & type=="indel" ~ "PASS",
            TRUE ~ "FILTERED")) %>%
        filter(filter=="PASS") %>%
        return
}

combine_datasets = function(ENST_CLASS_META_COUNT1,ENST_CLASS_META_COUNT2,TYPE="count"){
    if(TYPE=="count"){
        full_join(ENST_CLASS_META_COUNT1,ENST_CLASS_META_COUNT2,by=c("enst","class")) %>% 
            mutate_all(funs(ifelse(is.na(.),0,.))) %>% 
            mutate(
                cohort.wide.enst.wide.count=cohort.wide.enst.wide.count.x + cohort.wide.enst.wide.count.y,
                meta=map2(meta.x,meta.y,~rbind(.x,.y))) %>%
            select(c(enst,class,meta,cohort.wide.enst.wide.count)) %>% return
    }else if(TYPE=="rate"){
        full_join(ENST_CLASS_META_COUNT1,ENST_CLASS_META_COUNT2,by=c("enst","class")) %>% 
            mutate_all(funs(ifelse(is.na(.),0,.))) %>% 
            mutate( cohort.wide.enst.wide.rate=cohort.wide.enst.wide.rate.x + cohort.wide.enst.wide.rate.y ) %>%
            select(c(enst,class,cohort.wide.enst.wide.rate)) %>% return
    }
}

test_enrichment = function(YCU.COUNT,YCU.RATE,CLASSES.LIST,ENST__OF.INTEREST,GENOME.WIDE=FALSE){
    full_join(YCU.COUNT, YCU.RATE, by=c("enst","class")) %>%
        rename(exp=cohort.wide.enst.wide.rate,obs=cohort.wide.enst.wide.count) %>%
        mutate(obs=if_else(is.na(obs),0,as.numeric(obs))) %>% 
        mutate(exp=if_else(is.na(exp),0,as.numeric(exp))) %>% 
        filter(!is.na(exp)) %>%      
        nest(-enst) %>% 
        inner_join(ENST__OF.INTEREST) -> TMP

    TMP[rep(1:nrow(TMP), times=length(CLASSES.LIST)), ] %>% 
        mutate(classes=rep(CLASSES.LIST,each=nrow(TMP))) %>%
        mutate(data=map(data,make_LOF_row__one_ENST)) %>%
        mutate(data2=map2(data,classes,~summarise_classes(.x,.y))) %>%
        unnest(data2) %>% 
        select(-c(data,classes)) %>%
        left_join(enst_ddg2p) %>%
        { if(GENOME.WIDE){ group_by(.,class) %>% summarise(obs=sum(obs),exp=sum(exp)) %>% ungroup }else{ . } } %>%
        mutate(
            poisson.res=map(exp,~poisson.test(round(.x))),
            ci.low=map(poisson.res,~.x$conf.int[1]),
            ci.up=map(poisson.res,~.x$conf.int[2]),
            p.value=map2(obs,exp,~1-ppois(.x-1,.y))) %>%
        unnest(ci.low,ci.up,p.value) %>% return
}

summarise_classes = function(DT,CLASSES){
    CLASSES %>% unlist -> CLASSES.VEC
    tibble(
        class=paste(CLASSES.VEC,collapse="+"),
        meta=c(filter(DT, class %in% CLASSES.VEC) %>% select(meta)),
        obs=c(filter(DT, class %in% CLASSES.VEC) %>% pull(obs) %>% sum),
        exp=c(filter(DT, class %in% CLASSES.VEC) %>% pull(exp) %>% sum)
    ) %>% return
}

make_LOF_row__one_ENST = function(DT){
    bind_rows(DT,tibble(
            class=c(
                "LOF"
                ,"LOF_MIS"
                ,"LOF_CADDMIS"
                ,"CNV+LOF"
                ,"CNV+LOF_MIS"
                ,"CNV+LOF_CADDMIS"
            ),
            obs=c(
                filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS"))                   %>% pull(obs) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS","MIS"))             %>% pull(obs) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS","CADD.MIS"))         %>% pull(obs) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","CNV","NON","SPS","FS"))                   %>% pull(obs) %>% sum 
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","CNV","NON","SPS","FS","MIS"))             %>% pull(obs) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","CNV","NON","SPS","FS","CADD.MIS"))         %>% pull(obs) %>% sum
            ), 
            exp=c(
                filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS"))                   %>% pull(exp) %>% sum 
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS","MIS"))             %>% pull(exp) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","NON","SPS","FS","CADD.MIS"))         %>% pull(exp) %>% sum
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","CNV","NON","SPS","FS"))                   %>% pull(exp) %>% sum 
                ,filter(DT,class %in% c("SPAI.MIS","SPAI.SYN","CNV","NON","SPS","FS","MIS"))             %>% pull(exp) %>% sum
            )
        )
    )
}

plot_forrest = function(DF,CLASSES.OF.INTEREST,TITLE="aaa"){
    ggplot(DF,aes(x=class,y=exp)) +
        geom_point(size=1) +
        geom_point(aes(y=obs),size=1,color="red") +
        scale_x_discrete(limits=CLASSES.OF.INTEREST) +
        geom_errorbar(aes(ymax=ci.up,ymin=ci.low),width=0) +
        ggtitle(TITLE) 
}
