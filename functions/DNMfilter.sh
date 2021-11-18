# 1st arg: triodenovo result (Sample,chr,pos,ref,alt,whatever,,,,
# 2nd arg: ped
# 3rd arg: bam list (format: tab-delimited: SampleID PATH
# output: SampleID,chr,pos,score

# make DNM candidate file (comma-sep; numeric-sort; Family,chr(no "chr"header),pos
awk '{print $1","$2","$3}' $1|sed -e 's/Sample_//g' -e 's/chr//g' -e 's/X/999X/g' -e 's/Y/9999Y/g'  > tmp`basename $0`0.csv
grep -v -e MT -e GL -e hs tmp`basename $0`0.csv|sort -t"," -k1,1n -k2,2n -k3,3n|sed -e  's/999X/X/g' -e 's/9999Y/Y/g' > tmp`basename $0`01.csv
less tmpDNMfilter.sh0.csv| grep -e "MT"| sort -t"," -k1,1n -k2,2n -k3,3n > tmpDNMfilter.sh02.csv
less tmpDNMfilter.sh0.csv| grep -e "hs"| sort -t"," -k1,1n -k2,2n -k3,3n > tmpDNMfilter.sh04.csv
cat tmp`basename $0`01.csv tmp`basename $0`02.csv tmp`basename $0`04.csv|awk '{print "Sample_"$0}' > tmp`basename $0`.csv

# numeric-sort ped file
paste <(cut -f1,2 $2|sed 's/Sample_//g') $2|sort -k1,1n -k2,2n|cut -f3- > tmp`basename $0`.ped

# numeric-sort bam list
paste <(cut -f1 $3|sed 's/Sample_//') $3|sort -k1,1n|cut -f2- > tmp`basename $0`.list

# analysis
java -jar DNMFilter-0.1.1/DNMFilter.jar gbm \
    --reference hs37d5.fa \
    --pedigree tmp`basename $0`.ped \
    --bam tmp`basename $0`.list \
    --training DNMFilter-0.1.1/Training.264epi4ktrios.csv \
    --candidate tmp`basename $0`.csv \
    --configuration DNMFilter-0.1.1/Features.conf \
    --cutoff 0 \
    --output tmp`basename $0`3.csv > log`basename $0`.txt 2>&1 

# csv -> txt
awk -F"," 'BEGIN{OFS="\t"}{print $1,$2,$3,"NA","NA",$4}' tmp`basename $0`3.csv

