# 1st arg: vcf
# 2nd arg: ped file (tab-del: familyID, personalID, fatherID(or 0), motherID(or 0), sex(1:fa,2:mo), status(1:unaf,2:af))
# OUTPUT: score format (comma-del: sampleID,chr,pos,score)

TRIODENOVO="Triodenovo-0.05/bin/triodenovo"

# triodenovo analysis (option: --minDepth 0)
awk -F"\t" '$1~/^#/ || $9~/PL/ || $9~/GL/' <(less $1) > tmp`basename $0`.0.vcf 
$TRIODENOVO --ped $2  --minDQ -100 --minDepth 0 --in_vcf  tmp`basename $0`.0.vcf  --out_vcf tmp`basename $0`.vcf > log`basename $0`.txt 2>&1

# make list of proband
awk '$3!=0{print $2}' $2 > tmp`basename $0`Pt.txt

# vcf -> score format 
TriodenovoToScore.sh tmp`basename $0`.vcf tmp`basename $0`Pt.txt

