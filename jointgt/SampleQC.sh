VCF="vqsr_filter/cohort.recal.filtered.vcf"
SNPEFF="java -Xmx16g -jar snpEff.jar"
BCFTOOLS="bcftools"
GATK="java -Xmx64g -jar GenomeAnalysisTK.jar"
VCFTOOLS="vcftools"
REF="human_g1k_v37_fix.fasta"
RELATION="sample.relation.20191111.txt"
PED="trio_20191021.txt"
CAT="java  -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants"
SAMPLEYCU_POP=sample_country20191112.edit.txt
SAMPLE_PARENT1_PARENT2=sample_parent1_parent2.20191120.txt
SAMPLEHEALTHY=sample.healthy.invcf.txt
SAMPLEHEALTHYJPN=sample_phenotype__unrelated.jpn.healthy.txt
SAMPLEHEALTHYJPN200=sample_phenotype__unrelated.jpn.healthy200.txt
SAMPLEHEALTHYFOREIGN=sample.healthy.foreign.txt 
EXOMESUMMARY=exome_result.refGene.hg19.20200120_231439.txt
PC=/archive3/hamanaka/resource/variant.sanger.confirmed20200114.txt
export PERL5LIB=$PERL5LIB:/usr/local/genome/vcftools-0.1.17/share/perl5/
export LD_LIBRARY_PATH=/usr/local/genome/zlib-1.2.11/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/genome/samtools-1.9/bin/:$PATH
source Variable.sh
source VcfFunctions.sh

mv cohort.raw.vcf cohort.raw.original.vcf
$BCFTOOLS norm -m -any -f $REF cohort.raw.original.vcf > cohort.raw.norm.vcf 
make -f Makefile.vqsr.edit build vqsr > log.vqsr.txt 2>&1
cat <(head -n6 cohort.snp.recal.filtered.vcf)   <(head -n7 cohort.indel.recal.filtered.vcf|tail -n3) <(tail -n+7 cohort.snp.recal.filtered.vcf)   > cohort.snp.recal.filtered.edit.vcf
cat <(head -n7 cohort.indel.recal.filtered.vcf) <(head -n6 cohort.snp.recal.filtered.vcf  |tail -n2) <(tail -n+8 cohort.indel.recal.filtered.vcf) > cohort.indel.recal.filtered.edit.vcf
$CAT -R $REF -V cohort.snp.recal.filtered.edit.vcf -V cohort.indel.recal.filtered.edit.vcf -out cohort.recal.filtered.cat.vcf -assumeSorted > log.cat.txt
awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' cohort.recal.filtered.cat.vcf > n.vcf
cut -f1-10 n.vcf > s.vcf
for CHR in `seq 1 22` X Y; do $VCFTOOLS --vcf n.vcf --chr $CHR --out n.chr$CHR.vcf --recode --recode-INFO-all; done
awk '$1 ~ /^#/ || ($1 !~ /^#/ && $5!="*"){print $0}' cohort.recal.filtered.cat.site.vcf > s.vcf
make -f Makefile.edit Annovar > log.makefile.txt 2>&1 

$SNPEFF -v -onlyProtein -canon GRCh37.75 vqsr_filter/s.vcf                  > vqsr_filter/ss.vcf 2> log.snpeff.txt 
normalize_multi_anno vqsr_filter/ss.vcf |sort -k1,1n -k2,2n -u              > vqsr_filter/ssn.vcf
check_variant_sample_combi $PC n.vcf                                        > pc_status.txt
extract_genotype_from_variant_sample_combi $PC2 n.vcf                        > pc_gt.txt
extract_vqslod2 vqsr_filter/s.vcf                                           > id_vqslod.txt
extract_inbre vqsr_filter/s.vcf                                             > id_inbre.txt
extract_called_variant_sample_genotype_comb vqsr_filter/n.vcf 10            > id_sample_gt__ltAC10.txt 
check_vqslod_by_roc.R --pc pc_status.txt --vqslod id_vqslod.txt             --out vqslod.bin_sensitivity__snv.or.indel.png 
check_vqslod_by_roc.R --pc pc_status.txt --vqslod id_inbre.txt              --out inbre.bin_sensitivity__snv.or.indel.png 
extract_varid_ano_from_exomsummary $EXOMESUMMARY                            > id_annovar.txt
extract_varid_ano_from_snpeffvcf vqsr_filter/ss.vcf                         > id_snpeff.txt  
extract_varid_ano_from_snpeffvcf2 vqsr_filter/ss.vcf                        > id_snpeff2.txt 
awk -F"[\t|]" '{split($8,A,":"); print $1"-"$2"-"$3"-"$4"\t"$5"\t"$7"\t"A[1]}' ssn.vcf|sort -u > id_snpeff__normalized.txt 
select_sample_and_calc_af vqsr_filter/n.vcf.gz $SAMPLEHEALTHY               > id_maf.ycu.ctrl.txt
select_sample_and_calc_af vqsr_filter/n.vcf.gz $SAMPLEHEALTHYJPN            > id_maf.ycu.ctrl.jpn.txt
select_sample_and_calc_af vqsr_filter/n.vcf.gz $SAMPLEHEALTHYJPN200         > id_maf.ycu.ctrl.jpn.200.txt
select_sample_and_calc_af vqsr_filter/n.vcf.gz $SAMPLEHEALTHYFOREIGN        > id_maf.ycu.ctrl.foreign.txt
awk '$2<0.0001{print $1}' id_maf.ycu.ctrl.txt                               > id.rare.ycu.txt 
extract_rarevariantid_from_exomsummary $EXOMESUMMARY 0.00005                > id.rare.exac.txt

# sample qc
PaternityPihatMin=0.4
PaternityPihatMax=0.7
PaternityZ0Max=0.1
PihatRelatedness=0.125
ContamiPihatmeanThreshold=0.075
FMinFemale=-0.48 
FMaxFemale=0.48 
YMaxFemale=75
YMinMale=300 
FMinMale=0.75 
TiTvRatioMin=3.04
TiTvRatioMax=3.57
NSnpMin=4350 
NSnpMax=4900 
NInsMin=33 
NInsMax=60 
NDelMin=57
NDelMax=97 
InsDelRatioMin=0.4 
InsDelRatioMax=0.88 
HetHomRatioMin=1.8
HetHomRatioMax=5.2 
NPcPopUmap=7
PostJointGTIbd.sh $VCF $PaternityPihatMin $PaternityPihatMin $PaternityZ0Max             > sample1_sample2_ibd.data_ibd.status.txt
PostIbdRelatedness.sh $PihatRelatedness pair_z0_pihat.txt                                > sample_relatedsamples.txt
PostIbdContami.sh ../ibd/SMI contami.png $ContamiPihatmeanThreshold                      > sample_contami.data_contami.status.txt
PostJointGTSex.sh $VCF sex.png $FMinFemale $FMaxFemale $YMaxFemale $YMinMale $FMinMale   > sample_sex.data_sex.status.txt
PostJointGTPopulation.sh $VCF $SAMPLEYCU_POP $NPcPopUmap pop.umap.png                      sample_umaps_pop.rf.txt 
PostJointGTMetrics.sh $VCF titvratio.png count.png $TiTvRatioMin $TiTvRatioMax $NSnpMin $NSnpMax $NInsMin $NInsMax $NDelMin $NDelMax $InsDelRatioMin $InsDelRatioMax $HetHomRatioMin $HetHomRatioMax sample_titvratio.data_titvratio.status.txt sample_count.data_count.status.txt

merge_sampleqc.R \
    coerce/sample.invcf.txt \
    $SAMPLE_PARENT1_PARENT2 \
    coerce/sample_coerce.status.txt \
    ibd/sample1_sample2_ibd.data_ibd.status.txt \
    ibd/sample_relatedsamples.txt \
    ibd/sample_contami.data_contami.status.txt \
    sex/sample_sex.data_sex.status.txt \
    metrics/sample_titvratio.data_titvratio.status.txt \
    metrics/sample_count.data_count.status.txt \
    pop/sample_umaps_pop.rf.edit.txt \
    capture/sample_capture.umaps_capture.status.txt \
    sampleqc.20191202.txt

check_sampleqc_by_plot.R sampleqc.20191202.txt sampleqc.20191121.png
awk -F"\t" '$4=="ok" && $7=="ok" && $10!="outlier" && $19=="ok" && $21=="ok" && $23=="ok" && $25=="ok" && $27=="ok" && $29=="ok" && $42 !~ /unk/{print $1}' sampleqc.20191202.txt > sample.qced.txt
