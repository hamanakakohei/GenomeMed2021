# 1st arg: denovovariantcandidateid_pt_fa_mo.txt
# 2nd arg: vcf
# 3rd arg: id_snpeff.txt (id anno gene enst
# 4rd arg: ped
# 5th arg: output header
# output: DenovoFilter format for snv, indel, family 


awk 'NR==FNR{
        VAR_SAMP[$1";"$2]
        VAR_SAMP[$1";"$3]
        VAR_SAMP[$1";"$4]
    }NR!=FNR && $1=="#CHROM"{
        for(i=10;i<=NF;i++){SAMP[i]=$i}
    }NR!=FNR && $1 !~ /^#/{
        N=split($9,FORMAT,":")
        for(j=1;j<=N;j++){
            if(FORMAT[j]=="AD"){AdPos=j;break}
            if(j==N){next}
        }
        for(i=10;i<=NF;i++){
            VarThisRow=$1"-"$2"-"$4"-"$5
            SampThisCol=SAMP[i]
            split($i,Gt,":"); Ad=Gt[AdPos]
            if(VarThisRow";"SampThisCol in VAR_SAMP){
                print VarThisRow,SampThisCol,Ad
            }
        }
    }' $1 $2 > tmp`basename $0`.id_sample_depth.txt

denovofilterFormat2.R $1 tmp`basename $0`.id_sample_depth.txt $3 ${5}.snv.txt ${5}.indel.txt

awk '$3!=0{print $1,$2,$5}' $4|awk 'BEGIN{OFS="\t";print "family_id","individual_id","sex"}$3==1{$3="M"; print $0}$3==2{$3="F"; print $0}' > ${5}.family.txt
