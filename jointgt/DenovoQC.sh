VCFTOOLS="vcftools"
VCF="../vqsr_filter/n.vcf"
PED="trio_20191021.sexedit.txt"
PED2="trio.invcf.sex.ped"
SAMPLE_BAMPATH="trio_bampath.txt"
IDRAREEXAC=../id.rare.exac.txt
IDRAREYCU=../id.rare.ycu.txt
VCFRAREDNV=../vqsr_filter/nrd.vcf
SANGER=denovo.sanger.confirmed20191201.txt
OUTVARIANTQC=denovoqc.20200206.txt
OUTPNG=denovoqc.20200206.png
OUTSAMPLEQC=sampleqc.denovo.20200206.txt
export PERL5LIB=$PERL5LIB:/usr/local/genome/vcftools-0.1.17/share/perl5/
export LD_LIBRARY_PATH=/usr/local/genome/zlib-1.2.11/lib:$LD_LIBRARY_PATH
source VcfFunctions.sh

# trio sample in vcf?
join <(awk '{print $2,$1}' $PED|sort) <(head -n1000 $VCF|awk '$1 ~ /^#CHROM/'|cut -f10-|tr '\t' '\n'|sort)|cut -d" " -f2|sort|uniq -c|awk '$1==3{print $2}'|sort|join -t$'\t' <(sort $PED) - > trio.invcf.ped
join -1 2 -2 1 <(sort -k2,2 2536trio.invcf.ped) <(awk 'NR>1{print $1,$10}' ../sampleqc.20191121.txt|sed -e 's/female/2/' -e's/male/1/'|sort)|awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$4,$6}'|sort|awk '$3==0{print $0"\t"1}$3!=0{print $0"\t"2}' > trio.invcf.sex.ped

$VCFTOOLS --vcf $VCF --mendel $PED2 --out id.fam > log.mendel.txt 2>&1
awk '($6=="0\/1" || $6=="1\/0" || $6=="1\/1") && $7=="0\/0" && $8=="0\/0"{print $1"-"$2"-"$3"-"$4"\t"$5}' id.fam.mendel|\
    awk -F"_" '{print $1"_"$2"\t"$3"_"$4"\t"$5"_"$6}'|awk '$1 !~ /*/'                       > id.dnv_pt_fa_mo.txt
join <(sort $IDRAREEXAC) <(sort $IDRAREYCU)|sort|join -t$'\t' - <(sort id.dnv_pt_fa_mo.txt) > id.rare.dnv_pt_fa_mo.txt
cut -f1 id.rare.dnv_pt_fa_mo.txt                                                            > id.rare.dnv.txt
subset_vcf $VCF id.rare.dnv.txt                                                             > ../vqsr_filter/nrd.vcf
extract_vqslod2 $VCFRAREDNV                                                                 > id.rare.dnv_vq.txt 
TrioDenovo.sh $VCFRAREDNV $PED2                                                             > id.rare.dnv_td.txt
DNMfilter.sh id.rare.dnv_td.txt $PED2 $SAMPLE_BAMPATH                                       > id.rare.dnv_dn.txt
denovogear.sh $VCFRAREDNV $PED2                                                               id.rare.dnv_dg.txt
denovofilterFormat2.sh id.rare.dnv_pt_fa_mo.txt $VCFRAREDNV ../id_snpeff.txt $PED2 denovofilterFormat 
denovofilterFormat3.R denovofilterFormat.snv2.txt denovofilterFormat.indel2.txt               id.rare.dnv_df.txt 
merge_denovoqc.R \
    id.rare.dnv_pt_fa_mo.txt \
    id.rare.dnv_vq.txt \
    id.rare.dnv_td.txt \
    id.rare.dnv_dn.txt \
    id.rare.dnv_dg.txt \
    id.rare.dnv_df.txt \
    $SANGER \
    ../id_snpeff2.txt \
    ../id_annovar.txt \
    $OUTVARIANTQC 

check_denovoqc_by_plot.R \
    --denovoqc $OUTVARIANTQC \
    --ped $PED2 \
    --sampleqc ../sampleqc.20191202.txt \
    --snv_vq -8.18 \
    --snv_td 5.2 \
    --snv_dn 0.199 \
    --snv_dg 0.0025 \
    --indel_vq -3.19 \
    --indel_td 5.5 \
    --indel_df \
    --filteredtotal_max 10 \
    --out_png $OUTPNG \
    --out_txt $OUTSAMPLEQC 
