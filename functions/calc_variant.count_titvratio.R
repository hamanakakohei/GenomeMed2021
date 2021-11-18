#! /usr/bin/Rscript
titv_path=commandArgs(trailingOnly=TRUE)[1]
output.table_titv=commandArgs(trailingOnly=TRUE)[2]
output.plot_titv=commandArgs(trailingOnly=TRUE)[3]
count_path=commandArgs(trailingOnly=TRUE)[4]
output.table_count=commandArgs(trailingOnly=TRUE)[5]
output.plot_count=commandArgs(trailingOnly=TRUE)[6]
histogram_titvratio(titv_path,output.table_titv,output.plot_titv)
histogram_n.snv_n.ins_n.del_del.ins.ratio_het.homo.ratio(count_path,output.table_count,output.plot_count)

