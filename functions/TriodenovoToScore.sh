#Usage:
# 1st arg: Triodenovo-annotated vcf 
# 2nd arg: proband list 
# output: Sample,chr,pos,ref,alt,score

# Sample,chr,pos,ref,alt,score
NROW=`awk '$1=="#CHROM"{print NR}' $1`
awk '
    NR=='"${NROW}"' {
    for (i=10; i<=NF; i++) {
        SAMPLE[i]=$i
    }
}   NR>'"${NROW}"' {
    for (i=10; i<=NF; i++) {
        if ( $i !~ /^\./ ){
            split($i,ARRAY,":")
            print SAMPLE[i], $1, $2, $4, $5, ARRAY[2]
        }
    }
}' $1 > tmp`basename $0`.txt 

# only for proband
join <(sort tmp`basename $0`.txt) <(sort $2)|tr ' ' '\t' 

