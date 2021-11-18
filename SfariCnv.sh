DUPBED=/archive3/hamanaka/resource/SSC_SV/phase1/Project_REI_10816_B01_GRM_WGS.B38.StructuralVariants.2018-11-12/union_2.0.1696_batch21/REI_10816_B01_GRM_WGS_2.0.1696_batch21.05-07-2018.union.meta.del.merged.0.8.anno.genes.gFeatures.bed
DUPVCF=/archive3/hamanaka/resource/SSC_SV/phase1/Project_REI_10816_B01_GRM_WGS.B38.StructuralVariants.2018-11-12/union_2.0.1696_batch21/GS.sv.genotype.filtered.vcf
SV_PATH=/archive3/hamanaka/resource/SSC_SV/
FAM_STATUS_SAMPLE=/betelgeuse01/analysis/hamanaka/png/globus/nygc_sfari_id_map.csv
LIFTOVER="java -jar /usr/local/bio/src/picard-tools-2.10.0/picard.jar LiftOverIntervalList"
HG=/archive3/hamanaka/resource/hg19.fa.gz
DICT=/archive3/hamanaka/resource/hg19.fa.dict
CHAIN=/archive3/hamanaka/resource/hg38ToHg19.over.chain.gz
MAKEDICT="java -jar /usr/local/bio/src/picard-tools-2.10.0/picard.jar CreateSequenceDictionary"
export PATH=/betelgeuse01/analysis/hamanaka/function/:$PATH
source Variable.sh
source VcfFunctions.sh
export -f paste_sample_name_at_last
export -f spread_each_sample_calls
export -f output_format_types
export -f check_vcf_samples

find $SV_PATH -name "*.genomestrip.del.bed"|xargs -I{} sh -c "paste_sample_name_at_last {}" > merge.eachsamplebed.txt
find $SV_PATH -name "GS.sv.genotype.filtered.vcf"|grep -v "$DUPVCF"|xargs -I{} sh -c "spread_each_sample_calls {}" > merge.eachvcf.txt
find $SV_PATH -name "GS.sv.genotype.filtered.vcf"|xargs -I{} sh -c "output_format_types {}" > all.format.types.txt
find $SV_PATH -name "GS.sv.genotype.filtered.vcf"|grep -v "$DUPVCF"|xargs -I{} sh -c "check_vcf_samples {}" > path.vcf_sample.txt
find $SV_PATH -name "*gFeatures.bed"| grep -v "$DUPLICATED"| xargs -I{} cat {}|sed 's/^/chr/'|sort -u > all.batch.gFeatures.bed
spread_gfeatures_bed all.batch.gFeatures.bed > spread.jointbed.txt

CNV=cnv.id_count_gt.members_denovo.status__permissive.thr10000.txt
$MAKEDICT REFERENCE=$HG OUTPUT=$DICT
cat $DICT <(awk -F"[_\t]" 'BEGIN{OFS="\t"}NR>1 && !a[$1"-"$2"-"$3]++{print $1,$2,$3,"+",NR}' $CNV) > cnv.hg38.list
$LIFTOVER I=cnv.hg38.list O=cnv.hg19.list SD=$DICT CHAIN=$CHAIN > log.liftover.txt 2>&1
join -1 5 -2 5 <(awk '$1 !~ /^@/' cnv.hg38.list|sort -k5,5) <(awk '$1 !~ /^@/' cnv.hg19.list|sort -k5,5)|awk '{print $2"_"$3"_"$4,$6"_"$7"_"$8}'|sort|\
    join - <(tail -n+2 $CNV|sort)|tr ' ' '\t'|cut -f2-|cat <(head -n1 $CNV) - > cnv.id_count_gt.members_denovo.status.hg19.permissive.thr10000.txt
