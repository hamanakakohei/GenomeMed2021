GATK="java -jar GenomeAnalysisTK.jar"
REF=hs37d5.fa
VCFYCU="../vqsr_filter/cohort.recal.filtered.vcf"
PLINK="plink"
VCF1000G=1000g/ALL.autosome3.vcf
SAMPLE1000G_POP=1000g/1000genomes_20130606_sample_info.txt.gz
SAMPLEYCU_POP=sample_country20191112.txt

# select snp position in 1000 genome vcf
$GATK -T SelectVariants -R $REF -V $VCF1000G -conc $VCFYCU -o 1000gp.vcf > log1.txt 2>&1

# merge my vcf and 1000 genome vcf
merge_two_vcfs_for_common_variants $VCFYCU $VCF1000G > merge.vcf

#filter by MAF HWE etc
$PLINK --vcf merge.vcf --out merge.mafgenohwebi.vcf --mind 0.1 --maf 0.01 --geno 0.01 --hwe 0.05 --double-id --biallelic-only --recode > log4.txt 2>&1

#modify snpID
cp merge.mafgenohwebi.vcf.map merge.mafgenohwebi.vcf.map.tmp
awk '{OFS="\t"}{print $1,$1"_"$4,$3,$4}' merge.mafgenohwebi.vcf.map.tmp > merge.mafgenohwebi.vcf.map 

#indep
$PLINK --file merge.mafgenohwebi.vcf --out merge.mafgenohwebi.indep --indep 50 5 2 > log5.txt 2>&1

#indep extract
$PLINK --file merge.mafgenohwebi.vcf --recode transpose --out merge.mafgenohwebi.indep --extract merge.mafgenohwebi.indep.prune.in > log6.txt 2>&1

#genome
$PLINK --tfile merge.mafgenohwebi.indep --genome --out merge.mafgenohwebi.indep > log7.txt 2>&1

#PCA
$PLINK --tfile merge.mafgenohwebi.indep --out merge.mafgenohwebi.indep.pca \
    --read-genome merge.mafgenohwebi.indep.genome --cluster --pca 30  header > log9.txt 2>&1

# plot pc
barplot.R merge.mafgenohwebi.indep.pca.eigenval 30 merge.mafgenohwebi.indep.pca.eigenval.png

# umap
cat <(less $SAMPLE1000G_POP) <(tail -n+2 $SAMPLEYCU_POP) > sample_pop.txt
python3 umap.py merge.mafgenohwebi.indep.pca.eigenvec sample_pop.txt 7 sample_umaps_pop.txt umap.png
