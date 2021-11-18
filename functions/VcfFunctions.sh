extract_vqslod2(){ #1st arg: vcf; output: status"\t"vqslod score
    awk -F"\t" '$1 !~ /^#/{
        printf $1"\t"$2"\t"$4"\t"$5"\t"
        n=split($8,A,";")
        for(i=1;i<=n;i++){
            if(A[i] ~ /^VQSLOD=/){print A[i]; break}
            if(i == n){print "NoScore"}
        }
    }' <(less $1) | sed -e 's/VQSLOD=//g'
}

merge_two_vcfs_for_common_variants(){ #$1: 1st vcf; $2: 2nd vcf
    export PATH=/usr/local/genome/samtools-1.6/bin/:$PATH
    export PERL5LIB=$PERL5LIB:/usr/local/genome/vcftools-0.1.15/share/perl5/
    SAMTOOLS=/usr/local/genome/samtools-1.6/bin/
    bgzip tmp${FUNCNAME[0]}4.vcf
    tabix tmp${FUNCNAME[0]}4.vcf.gz
    bgzip tmp${FUNCNAME[0]}5.vcf
    tabix tmp${FUNCNAME[0]}5.vcf.gz
    vcf-merge tmp${FUNCNAME[0]}4.vcf.gz tmp${FUNCNAME[0]}5.vcf.gz 
}


