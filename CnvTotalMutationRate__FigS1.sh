BEDTOOLS="bedtools intersect"
SV="gnomad_v2_sv.sites.bed.gz"
PLOF="../exon_for_plof.bed"
CG="../transcript_for_cg.bed"
ENST_MARGIN="../enst_margin.1000k.0309.bed"
ASD="dnSV_ASD.bed"

# annotate gnomad bed with plof or cg 
ENST_MARGIN="../enst_margin.1000k.0309.bed"
less $SV|awk 'NR>1 && $5=="DEL"'|$BEDTOOLS -a stdin -b $PLOF -wa -wb|awk -F"\t" '!a[$4"-"$NF]++'|$BEDTOOLS -a stdin -b $ENST_MARGIN -wa -wb -f 1|awk -F"\t" '($3-$2)<1000000 && $NF==$(NF-4)' > gnomad_plof.1000k.0309.bed
less $SV|awk 'NR>1 && $5=="DUP"'|$BEDTOOLS -a stdin -b $CG -wa -wb -F 1|$BEDTOOLS -a stdin -b $ENST_MARGIN -wa -wb -f 1|awk -F"\t" '($3-$2)<1000000 && $NF==$(NF-4)' > gnomad_cg.1000k.0309.bed
awk -F"\t" '{print $NF}' gnomad_plof.1000k.0309.bed|sort|uniq -c > count_plof_enst_mafall.1000k.0309.txt
awk -F"\t" '{print $NF}' gnomad_cg.1000k.0309.bed|sort|uniq -c > count_cg_enst_mafall.1000k.0309.txt
awk -F"\t" '{a[$4]=a[$4]$NF";"}END{for(i in a)print i"\t"a[i]}' gnomad_plof.1000k.0309.bed > sv.name_plof.enst.1000k.0309.txt
awk -F"\t" '{a[$4]=a[$4]$NF";"}END{for(i in a)print i"\t"a[i]}' gnomad_cg.1000k.0309.bed > sv.name_cg.enst.1000k.0309.txt

# annotate ASD paper (Werling et al
ENST_MARGIN="../enst_margin.1000k.0309.bed"
less $ASD|awk '$5=="DEL"'       |$BEDTOOLS -a stdin -b $PLOF -wa -wb|awk -F"\t" '!a[$4"-"$NF]++'|$BEDTOOLS -a stdin -b $ENST_MARGIN -wa -wb -f 1|awk -F"\t" '($3-$2)<1000000 && $NF==$(NF-4)' > asd_plof.bed
less $ASD|awk '$5=="DUP"'       |$BEDTOOLS -a stdin -b $CG -wa -wb -F 1|$BEDTOOLS -a stdin -b $ENST_MARGIN -wa -wb -f 1|awk -F"\t" '($3-$2)<1000000 && $NF==$(NF-4)' > asd_cg.bed
awk -F"\t" '$6 ~ /p/{print $NF}' asd_plof.bed|sort|uniq -c > count_plof_enst_asd.txt
awk -F"\t" '$6 ~ /p/{print $NF}' asd_cg.bed|sort|uniq -c > count_cg_enst_asd.txt
cat asd_plof.bed asd_cg.bed|awk -F"\t" '{a[$4"-"$5"-"$6]=a[$4"-"$5"-"$6]$NF";"}END{for(i in a)print i"\t"a[i]}'|tr '-' '\t' > asd_enst.txt
