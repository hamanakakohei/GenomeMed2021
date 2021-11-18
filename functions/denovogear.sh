#1st: vcf
#2nd: ped
#3rd: output

source /usr/local/genome/denovogear-1.1.1/env.sh
dng dnm auto --rd_cutoff 1 --ped $2 --vcf $1 --pp_cutoff 0.0000000000001|awk 'NF!=0 && NF!=3 && NF!=9' > tmp`basename $0`.txt

denovogear.R tmp`basename $0`.txt $3

