# input: vcf (automatically extract biallelic variant here
# output: table of depth (DP)

NROW=`head -n 1000 $1|awk '$1=="#CHROM"{print NR}'`

awk 'BEGIN{OFS="\t"} NR=='"${NROW}"' {
    printf "VariantID\t" 
    for (i=10; i<=(NF-1); i++) { printf $i"\t" }
                    
    ## for last column
    printf $NF"\n"
} NR>'"${NROW}"' && $4 !~ /,/ && $5 !~ /,/ && $9 ~ /DP/{    
    printf $1"-"$2"-"$4"-"$5"\t"
    NOFFORMAT=split($9,FORMAT,":")
        
    for(i=1;i<=NOFFORMAT;i++){
        if(FORMAT[i]=="DP"){
            for(j=10;j<=(NF-1);j++){
                split($j,GENOTYPE,":")
                printf GENOTYPE[i]"\t"
            }

            # for last col.
            split($NF,GENOTYPE,":")
            printf GENOTYPE[i]"\n"
            break
        }
    }
}' $1

