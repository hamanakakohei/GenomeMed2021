gff3_path="gencode.v19.annotation.gff3.gz"
canonical_path="gnomad.v2.1.1.canonical.txt"

less $gff3_path|tail -n+8|cut -f9 > gff3_an.txt

less $gff3_path|tail -n+8|cut -f3 > gff3_an_type.txt

less $gff3_path|tail -n+8|awk -F"\t" '{print $1"\t"$4"\t"$5}' > gff3_an_start_end.txt

cut -d";" -f4 gff3_an.txt|cut -d"=" -f2|cut -d"." -f1 > gff3_an_enst.txt

cut -d";" -f2 gff3_an.txt|cut -d"=" -f2|cut -d"." -f1 > gff3_an_gene.txt

awk -F";" '{
    for(i=1;i<=NF;i++){
        if(i==1){
            gene_type="NA";transcript_type="NA"
        }
        if($i ~ /gene_type/){gene_type=$i}
        if($i ~ /transcript_type/){transcript_type=$i}
        if(i==NF){
            print gene_type"\t"transcript_type
        }
    }   
}' gff3_an.txt|tr '\t' ';' > gff3_an_att.txt

paste gff3_an_type.txt gff3_an_enst.txt gff3_an_att.txt|grep -f AttOI.txt|awk '$1=="transcript"{print $2}'|sort -u > enst_protein.coding.txt 

join <(sort $canonical_path) <(sort enst_protein.coding.txt) > enst_canonical.protein.coding.txt

paste gff3_an_type.txt gff3_an_start_end.txt gff3_an_enst.txt|awk '$1=="transcript"'|cut -f2- > enst_start_end.txt
