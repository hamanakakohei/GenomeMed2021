PLINK=plink
VCF="../vqsr_filter/cohort.recal.filtered.vcf"
GATK="java -jar GenomeAnalysisTK.jar"
REF="hs37d5.fa"

$PLINK --vcf $VCF --out p --maf 0.1 --geno 0.1 --hwe 0.05 --biallelic-only --recode vcf --double-id > log.txt 2>&1
$GATK -T VariantEval -R $REF -eval p.vcf -o titvratio.txt -noST -noEV -EV TiTvVariantEvaluator -ST Sample > log.titvratio.txt 2>&1
$GATK -T VariantEval -R $REF -eval p.vcf -o count.txt -noST -noEV -EV CountVariants -ST Sample > log.count.txt 2>&1
awk 'BEGIN{OFS="\t"}{$1=$1; print $0}' titvratio.txt > titvratio2.txt
awk 'BEGIN{OFS="\t"}{$1=$1; print $0}' count.txt > count2.txt
calc_variant.count_titvratio.R titvratio2.txt sample_titvratio.txt titvratio.png count2.txt sample_nSNPs_nInsertions_nDeletions_insertionDeletionRatio_hetHomRatio.txt count.png

# status annotation
awk 'BEGIN{OFS="\t"}NR==1{print $0}NR>1{n=split($1,A,"_");if(n==2){print A[1],$2
    }else if(n==4){print A[1]"_"A[2],$2
    }else if(n==6){print A[1]"_"A[2]"_"A[3],$2}}' sample_titvratio.txt| \
        awk 'NR==1{print $0"\ttitvratio.status"}NR>1{if($2>3 && $2<3.6){print $0"\tok"}else{print $0"\toutlier"}}' > sample_titvratio.data_titvratio.status.txt

awk 'BEGIN{OFS="\t"}NR==1{print $0}NR>1{n=split($1,A,"_");if(n==2){print A[1],$2,$3,$4,$5,$6
    }else if(n==4){print A[1]"_"A[2],$2,$3,$4,$5,$6
    }else if(n==6){print A[1]"_"A[2]"_"A[3],$2,$3,$4,$5,$6}}' sample_nSNPs_nInsertions_nDeletions_insertionDeletionRatio_hetHomRatio.txt| \
        awk 'BEGIN{OFS="\t"}NR==1{print $1"\tnsnp\tnsnp.status\tnins\tnins.status\tndel\tndel.status\tinsdelratio\tinsdelratio.status\thethomratio\thethomratio.status"}NR>1{
            if($2>3750 && $2<5500){printf $1"\t"$2"\tok\t"}else{printf $1"\t"$2"\toutlier\t"};
            if($3>25 && $3<75){printf $3"\tok\t"}else{printf $3"\toutlier\t"};
            if($4>40 && $4<110){printf $4"\tok\t"}else{printf $4"\toutlier\t"};
            if($5>0.35 && $5<1){printf $5"\tok\t"}else{printf $5"\toutlier\t"};
            if($6>1.8 && $6<5.3){printf $6"\tok\n"}else{printf $6"\toutlier\n"};
        }' > sample_count.data_count.status.txt


