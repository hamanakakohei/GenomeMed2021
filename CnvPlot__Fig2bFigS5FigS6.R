library(tidyverse)
p_path              = "enrich/prep_for_qqman_ycu.txt"
hummut_plof_path    = "enrich/ycu_plof_small.bed" 
hummut_cg_path      = "enrich/ycu_cg_small.bed"     
gene.region_path    = "gencode.v19.annotation.gff3.gz"
chr.df_path         = "enrich/chromosome_range.txt"
loeuf_path          = "gnomad.v2.1.1.lof_metrics.by_transcript.txt"
zs.v4v5v6_path      = "v4v5v6.t.txt.gz"
sample.healthy_path = "sample.healthy.txt"
dncnv.region__path  = "dncnv__ycu.for.gviz3.txt"
sample__rm__path    = "sample__RemovedForGviz.txt"
p_thr               = 0.00001 
sample__rm	    = read_tsv(sample__rm__path)
sample.healthy      = read_tsv(sample.healthy_path,col_names="sample") %>% pull(sample)
zs.v4v5v6           = read_tsv(zs.v4v5v6_path) %>% separate(Matrix,into=c("chr","start","end")) %>% mutate(chr=paste0("chr",chr))
zs.v4v5v6.sample    = colnames(select(zs.v4v5v6,starts_with("Sample_")))
hm_plof             = read_tsv(hummut_plof_path,col_types="cddcd",col_names=c("chr","sta","end","pt","copy")) %>% distinct(chr,sta,end,pt,copy) %>% mutate(chr=paste0("chr",chr))
hm_cg               = read_tsv(hummut_cg_path,col_types="cddcd",col_names=c("chr","sta","end","pt","copy")) %>% distinct(chr,sta,end,pt,copy) %>% mutate(chr=paste0("chr",chr))
chr.df              = read_tsv(chr.df_path)
p                   = read_tsv(p_path) %>% arrange(chr,start)
p_lof               = filter(p,type=="PLOF")
p_cg                = filter(p,type=="CG")
loeuf               = read_tsv(loeuf_path) %>% dplyr::select(c(chromosome,start_position,end_position,oe_lof_upper,canonical,transcript)) %>% 
	                    filter(canonical=="TRUE") %>% nest(-transcript) %>% mutate(points=map(.$data,interval_to_points)) %>% unnest(points) %>% mutate(chr=paste0("chr",chr))
type_cnv.region = p %>% nest(-type) %>% 
    mutate(
    data_segment=map(data,~give_segment_count_vc(.,p_thr))
    ,segment_region=map(data_segment,~give_segment_region(.))
    ,cnv.region=map2(segment_region,as.list(type),segment_to_overlap.cnv.region)
    ,segment_region=map2(segment_region,cnv.region,~ cbind(.x,sample=.y$sample) %>% as_tibble %>% mutate(sample=as.character(sample)))
)
sample.random.v4v5v6 = sample(intersect(sample.healthy,zs.v4v5v6.sample),500)

library(Gviz)
gtrack <- GenomeAxisTrack()
grtrack = readRDS("enrich/grtrack.gencode.rds")
zs.v4v5v6.gr = makeGRangesFromDataFrame(zs.v4v5v6,keep.extra.columns=T) 
atrack_YCU_PLOF = AnnotationTrack(start=hm_plof$sta,end=hm_plof$end,chromosome=hm_plof$chr,group=hm_plof$pt,genome="hg19",fill="red",name="LOF")
atrack_YCU_CG = AnnotationTrack(start=hm_cg$sta,end=hm_cg$end,chromosome=hm_cg$chr,group=hm_cg$pt,genome="hg19",fill="blue",name="CG")
dtrack_P_PLOF = DataTrack(data=-log10(p_lof$pval),start=p_lof$start,end=p_lof$end,chromosome=p_lof$chr,genome="hg19",name="p-value")
dtrack_P_CG = DataTrack(data=-log10(p_cg$pval),start=p_cg$start,end=p_cg$end,chromosome=p_cg$chr,genome="hg19",name="-log10(p)")
options(ucscChromosomeNames=FALSE)
dtrack_LOEUF = DataTrack(data=loeuf$loeuf,start=loeuf$start,end=loeuf$end,chromosome=loeuf$chr,genome="hg19",name="LOEUF",ylim=c(0,1.7))
sapply(paste0("chr",c(1:22,"X","Y")),function(x)eval(parse(text=paste0("itrack_",x,"<<-readRDS('",x,".rds')"))))
map2(type_cnv.region$cnv.region,as.list(type_cnv.region$type),~plot_cnv_wrap(df=.x,type=.y,RATIO=30  ,NAME="_20210311ycu-1.png"))
map(list("PLOF","CG"),~plot_chr_wrap(chr.df=chr.df,type=.x))

# functions
segment_to_overlap.cnv.region = function(df.segment,type){
	library(GenomicRanges)
	if(type=="PLOF"){df.cnv.gr = GRanges(hm_plof$chr,IRanges(hm_plof$sta,hm_plof$end),sample=hm_plof$pt)
		}else if(type=="CG"){df.cnv.gr = GRanges(hm_cg$chr,IRanges(hm_cg$sta,hm_cg$end),sample=hm_cg$pt)}
	df.segment.gr = GRanges(df.segment$chr,IRanges(df.segment$start,df.segment$end,names=df.segment$segment.id))
	findOverlaps(df.segment.gr,df.cnv.gr) -> segment.id_hit.id
	data.frame(df.cnv.gr) %>% as_tibble %>% rownames_to_column(var="hit.id") %>% mutate(hit.id=as.integer(hit.id)) %>%
		right_join(tibble(segment.id=segment.id_hit.id@from,hit.id=segment.id_hit.id@to),by="hit.id") %>%
		nest(-segment.id) %>%
		mutate(cnv.region=map(.$data,function(df){
			mutate(df,start.cnv=min(start),end.cnv=max(end),sample=paste(sample,collapse=";")) %>% 
				dplyr::select(seqnames,start.cnv,end.cnv,sample) %>% 
				distinct %>% 
				return})) %>%
		unnest(cnv.region) %>% 
		dplyr::select(-data) %>% 
		dplyr::rename(chr=seqnames,start=start.cnv,end=end.cnv) %>% 
		dplyr::mutate(chr=as.character(chr)) %>%
        return
}	

plot_cnv_wrap = function(df,type,RATIO,NAME){sapply(1:nrow(df),plot_cnv_pt.separate,df,type,RATIO,NAME)}
plot_chr_wrap = function(chr.df,type){sapply(1:nrow(chr.df),plot_chr,chr.df,type)}

interval_to_points = function(df){
	points = seq(df$start_position,df$end_position,length=1000)
	tibble(chr=df$chromosome,start=points,end=points,loeuf=df$oe_lof_upper) %>% return
}

plot_cnv_pt.separate = function(segment.id,df,type,RATIO,NAME="_20210216ycu-6.png"){
    chr=df$chr[segment.id]
    start=df$start[segment.id]
    end=df$end[segment.id]
    samples.pt = df$sample[segment.id] %>% str_split(pattern=";") %>% unlist
    margin= (end-start)/RATIO
    eval(parse(text=paste0("itrack=itrack_",chr)))
    eval(parse(text=paste0("dtrack_P=dtrack_P_",type)))
    eval(parse(text=paste0("atrack_YCU=atrack_YCU_",type)))
    sapply(samples.pt,pt_to_correspo.kit.dtrack,zs.v4v5v6.gr,sample.random.v4v5v6,chr,start,end,margin) -> list.dtrack
    png(paste0(type,"segment",segment.id,NAME),width=9.45,height=9,units="cm",res=900)
    try(plotTracks(c(list(itrack,gtrack,grtrack,dtrack_LOEUF,atrack_YCU),list.dtrack) 
        ,chromosome=chr
    	,from=start-margin
    	,to=end+margin
        ,collapseTranscripts="longest"
	,just.group="above"
        ,sizes=c(0.5,1,4,1,0.5,1)
    ))
    dev.off()
}

pt_to_correspo.kit.dtrack = function(pt,zs.v4v5v6.gr,sample.random.v4v5v6,chr,start,end,margin){
    mcols(zs.v4v5v6.gr) %>% colnames -> sample.v4v5v6   
    if(pt %in% sample.v4v5v6){
        tmp.dtrack = make_dtrack.zs(zs.v4v5v6.gr,sample.random.v4v5v6,pt,chr,start,end,margin)
    }
    return(tmp.dtrack)
}

make_dtrack.zs = function(zs.gr,samples.random,samples.pt,chr.plot,start.plot,end.plot,margin.plot){
    zs.gr2 = zs.gr[seqnames(zs.gr)==chr.plot & start(zs.gr)>start.plot-margin.plot & end(zs.gr)<end.plot+margin.plot]
    zs.gr3 = zs.gr2[,c(samples.random,samples.pt)]
    DataTrack(zs.gr3,groups=1:(mcols(zs.gr3) %>% ncol),col=make_cols(zs.gr3,samples.pt),lwd=3,type="l",legend=FALSE) %>% return
}

make_cols = function(kit.gr,samples.pt){
	mcols(kit.gr) %>% colnames -> samples.kit
    as.character(samples.kit %in% samples.pt) %>% 
		str_replace(pattern="TRUE",replacement=rgb(1,0,0)) %>%
		str_replace(pattern="FALSE",replacement=rgb(0,0,0,alpha=0.05)) %>%
		return
}

plot_chr = function(chr.id,chr.df,type){
    chr=chr.df$chr[chr.id]
    start=chr.df$start[chr.id]
    end=chr.df$end[chr.id]
    eval(parse(text=paste0("itrack=itrack_",chr)))
    eval(parse(text=paste0("dtrack=dtrack_",type)))
    eval(parse(text=paste0("atrack=atrack_",type)))
    png(paste0(type,"chr",chr.id,"0808.png"))
    plotTracks(list(itrack,gtrack,dtrack,atrack)
        ,chromosome=chr
    	,from=start
    	,to=end
        ,ylim=c(0,20)
    )
    dev.off()
}

give_segment_region = function(df){
	vc.n=max(df$segment,na.rm=T)
	chr_vc=start_vc=end_vc=enst_vc=rep(NA,vc.n)
	for(i in 1:nrow(df)){
		id=df$segment[i]
		if(!is.na(id)){
			chr_vc[id]=df$chr[i]	
			start_vc[id]=min(start_vc[id],df$start[i],na.rm=T)
			end_vc[id]=max(end_vc[id],df$end[i],na.rm=T)	
			enst_vc[id]=paste0(enst_vc[id],";",df$enst[i])	
		}
	}
	tibble(segment.id=1:vc.n,chr=chr_vc,start=start_vc,end=end_vc,enst=enst_vc) %>% return
}

give_segment_count_vc = function(dt,p_thr){
	segment_count=0
	segment_vc=rep(NA,nrow(dt))
	segment_chr=dt$chr[1]
	for(i in 1:nrow(dt)){
		if(i==1){
			if(dt$pval[i]<p_thr){
				segment_count=segment_count+1
				segment_vc[i]=segment_count
			}
		}else{
			if(dt$pval[i]<p_thr){
				if(segment_chr==dt$chr[i]){
					if(is.na(segment_vc[i-1])){
						segment_count=segment_count+1
						segment_vc[i]=segment_count
					}else{
						segment_vc[i]=segment_count
					}	
				}else{
					segment_count=segment_count+1
					segment_vc[i]=segment_count
					segment_chr=dt$chr[i]
				}
			}
		}
	}
	cbind(dt,segment=segment_vc) %>% as_tibble %>% return
}

