# 1st arg: vcf 
# 2: PaternityPihatMin
# 3: PaternityPihatMax
# 4: PaternityZ0Max
# output: sample1 sample2 ibd.z0 ibd.pihat ibd.status(="paternal" or "outlier")  
VCFTOOLS="vcftools"
REF="hs37d5.fa"
PLINK="plink"
RELATION="sample.relation.20191111.txt"
VCF=$1
PaternityPihatMin=$2
PaternityPihatMax=$3
PaternityZ0Max=$4

#filter by MAF HWE etc
$PLINK --vcf-filter --vcf $VCF --out SM --mind 0.1 --maf 0.1 --geno 0.1 --hwe 0.05 --double-id --biallelic-only --recode > log2.txt 2>&1

#modifysnpID
cp SM.map SM.map.tmp
awk '{OFS="\t"}{print $1,$1"_"$4,$3,$4}' SM.map.tmp > SM.map 

#independent
$PLINK --file SM --out SMI --indep 50 5 2 > log3.txt 2>&1

#independent extract
$PLINK --file SM --recode transpose --out SMI --extract SMI.prune.in > log4.txt 2>&1

#genome
$PLINK --tfile SMI --genome --out SMI > log5.txt 2>&1

# plot the result
awk 'BEGIN{print "pair\tz0\tpihat"}NR>1{print $2";"$4"\t"$7"\t"$10}' SMI.genome > pair_z0_pihat.txt
awk 'BEGIN{print "pair\trelation"}{print $1";"$2"\t"$3"\n"$2";"$1"\t"$3}' $RELATION > pair_relation.txt
PlotIbd2.R pair_z0_pihat.txt pair_relation.txt ibd.relation.png > log6.txt 2>&1

# change for the last merge
awk -F"[\t;]" 'BEGIN{OFS="\t";print "sample1","sample2","ibd.z0","ibd.pihat"}NR>1{print $1,$2,$3,$4"\n"$2,$1,$3,$4}' pair_z0_pihat.txt|\
    awk -v pmin="${PaternityPihatMin}" -v pmax="${PaternityPihatMax}" -v zmax="${PaternityZ0Max}" 'NR==1{print $0"\tibd.status"}NR>1{
        if($4>pmin && $4<pmax && $3<zmax){print $0"\tpaternal"}else{print $0"\toutlier"}
    }' 

