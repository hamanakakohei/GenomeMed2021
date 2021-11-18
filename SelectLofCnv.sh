BEDTOOLS="bedtools intersect"
PLOF="../exon_for_plof.bed"
YCU="invcf_cnv.in.xcnv_n.target_sqs.txt"
TRIOQC="sampleqc.denovo.20191202.txt"

# ycu call
awk -F"\t" '$5>90 && $7>=4 && $8<90 && $9<90 && $10<=2{print $6,$1,$4}' $YCU|\
    awk -F"[:-]" '{print $1,$2,$3,$4,$5}'|\
    sed -e 's/DEL/1/g' -e 's/DUP/3/'|\
    sort -k4,4| join -1 4 -2 1 - <(awk 'NR>1 && $2=="ok"{print $1}' $TRIOQC|sort)| \
    awk 'BEGIN{OFS="\t";print "chr\tstart\tend\tpt\tcopy"}{print $2,$3,$4,$1,$5}' \
    > YcuCnv.bed
    
less YcuCnv.bed|awk 'NR>1 && $5<2 && ($3-$2)<1000000'  |$BEDTOOLS -a stdin -b $PLOF -wa -wb    |awk '!a[$1"-"$2"-"$3"-"$4"-"$5"-"$NF]++' > YcuCnv_plof_small.bed
less YcuCnv.bed|awk 'NR>1 && $5>2 && ($3-$2)<1000000'  |$BEDTOOLS -a stdin -b $CG -wa -wb -F 1 |awk '!a[$1"-"$2"-"$3"-"$4"-"$5"-"$NF]++' > YcuCnv_cg_small.bed
awk -F"\t" '{print $NF}' YcuCnv_plof_small.bed     |sort|uniq -c > count_plof_enst_YcuCnv_small.txt
awk -F"\t" '{print $NF}' YcuCnv_cg_small.bed       |sort|uniq -c > count_cg_enst_YcuCnv_small.txt

