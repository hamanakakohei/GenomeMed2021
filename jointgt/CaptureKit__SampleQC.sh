VCF="../vqsr_filter/cohort.snp.recal.filtered.vcf" 
NOFSAMPLE=13847
NOFPC=8

# vcf -> DP table
VcfToDepth.sh $VCF > DP.txt

TMP=`expr 1 + $NOFSAMPLE`
awk 'NF=='"${TMP}"'' DP.txt > DP2.txt
grep -v "\." DP2.txt > DP3.txt

# normalize depth with total depth 
DepthNormalization.sh DP3.txt > DP4.txt

# check very low depth sample
barplot2.R tmpDepthNormalization.sh.txt total.depth.png
awk 'BEGIN{OFS="\t"}NR==1{print "sample","depth","depth.status"}NR>1{if($2<10000000){print $1,$2,"lowdepth"}else{print $1,$2,"ok"}}' tmpDepthNormalization.sh.txt > sample_depth.data_depth.status.txt

# PC plot and result => 8, 13, or 17 PCs tsukau
PC_PLOT.R DP4.txt 30 capture.pc.png sample_capture.pcs.txt

awk -v NPC="${NOFPC}" 'NR==1{
    for(i=0;i<=NPC;i++){
        if(i==0){printf "sample\t"
        }else if(i==NPC){print "cap.pc"i
        }else{printf "cap.pc"i"\t"}
    }
}NR>1{
    for(i=0;i<=NPC;i++){
        if(i==NPC){j=i+1; print $j 
        }else{j=i+1; printf $j"\t"}
    }
}' sample_capture.pcs.txt > sample_cap.pcs.txt


# umap
SAMPLE_CAPTURE="sample_capture20191113.txt"
python3 umap2.py sample_capture.pcs.txt $SAMPLE_CAPTURE $NOFPC capture.batch_sample_umaps.txt capture.batch.umap.png

