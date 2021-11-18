# 1st arg: table of depth (header: Variant_ID, Sample_ID, Sample_ID,,,,
# 2nd arg: sample list you want to remove
# output: table of depth normalized with total depth

# calcu. total depth of each sample (1s col: Sample_ID, 2nd col: total depth
awk 'NR==1{
    for(i=2;i<=NF;i++){
        SAMPLE[i]=$i
    }
}NR>1{
    for(i=2;i<=NF;i++){
        TOTAL[SAMPLE[i]]+=$i
    }
}END{
    print "Sample_ID","Total_Depth"
    for(i=2;i<=NF;i++){
        print SAMPLE[i],TOTAL[SAMPLE[i]]
    }
}' $1 > tmp`basename $0`.txt

# normalize with total depth (same format as 1st arg
awk 'BEGIN{OFS="\t"}NR==FNR && FNR>1{
    TOTAL[NR]=$2
}NR!=FNR && FNR==1{
    print $0
}NR!=FNR && FNR>1{
    for(i=2;i<=NF;i++){
        $i=int($i*100000000/TOTAL[i])
    }
    print $0
}' tmp`basename $0`.txt $1 > tmp`basename $0`2.txt 

# randomly select and reduce the file size
awk 'NR % 20 == 1' tmp`basename $0`2.txt

