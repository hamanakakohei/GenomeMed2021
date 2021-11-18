gff3_path="gencode.v19.annotation.gff3.gz"

less $gff3_path|tail -n+8|cut -f9 > gff3_an.txt

less $gff3_path|tail -n+8|cut -f3 > gff3_an_type.txt

less $gff3_path|tail -n+8|awk -F"\t" '{print $1"\t"$4"\t"$5}' > gff3_an_start_end.txt

cut -d";" -f4 gff3_an.txt|cut -d"=" -f2|cut -d"." -f1 > gff3_an_enst.txt

cut -d";" -f2 gff3_an.txt|cut -d"=" -f2|cut -d"." -f1 > gff3_an_gene.txt

less $gff3_path|tail -n+8|paste gff3_an_enst.txt - |awk '$4=="exon"{LEN[$1]=LEN[$1]+$6-$5}END{for(i in LEN)print i,LEN[i]}' > enst_length.txt

awk -F";" '{
    for(i=1;i<=NF;i++){
        if(i==1){gene_type="NA";transcript_type="NA"}
        if($i ~ /gene_type/){gene_type=$i}
        if($i ~ /transcript_type/){transcript_type=$i}
        if(i==NF){print gene_type"\t"transcript_type}
    }   
}' gff3_an.txt|tr '\t' ';' > gff3_an_att.txt

paste gff3_an_type.txt gff3_an_gene.txt gff3_an_enst.txt gff3_an_att.txt|grep -f AttOI0130.txt|awk '$1=="transcript"{print $3,$2}'|sort|join - <(sort enst_length.txt)| \
    awk '{
        if(!($2 in ENST)){
            ENST[$2]=$1
            LEN[$2]=$3
        }else if(LEN[$2]<$3){
            ENST[$2]=$1
            LEN[$2]=$3
        }
    }END{
        for(i in ENST)print i,ENST[i],LEN[i]
    }'> ensg_enst_length0130.txt 

grep ENSGR ensg_enst_length0130.txt|awk '{print $1}' > ensgr
sed 's/ENSGR/ENSG0/g' < ensgr > ensg0
grep -wv -f ensgr -f ensg0 ensg_enst_length0130.txt > ensg_enst_length.nopar0130.txt
awk '{print $2}' ensg_enst_length.nopar0130.txt|sort -u > enst_nc0130.txt
join -v 1 <(sort enst_nc0130.txt) <(awk '$1=="chrX" || $1=="chrY"{print $4}' enst_start_end.txt|sort) > enst_nc0130.noXY.txt

paste gff3_an_type.txt gff3_an_start_end.txt gff3_an_enst.txt|awk '$1=="transcript"'|cut -f2-|awk '!a[$4]++' > enst_start_end.txt
