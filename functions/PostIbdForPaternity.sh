# 1st arg: ibd.genome file
# 2nd arg: ped file (no header; tab-del: FID, IID, FatherID (or 0), MotherID (or 0),,,,,, 
# 3rd arg: output name of plot (ex. Paternity.png


# from ibd file, make sample-sample (child-parent or parent-child) pair, Z0, PI_HAT 
awk '{print $1"-"$3,$7,$10}' $1 > tmp`basename $0`.txt

# from ped file, make child-parent pair
awk '$3!=0{print $2"-"$3"\n"$2"-"$4"\n"$3"-"$2"\n"$4"-"$2}' $2 > tmp`basename $0`2.txt

# make child-parent pair, Z0, PI_HAT file
join <(sort tmp`basename $0`.txt) <(sort tmp`basename $0`2.txt)|cat <(echo "pair Z0 PI_HAT") - > tmp`basename $0`3.txt

# make Z0, PI_HAT only file
cut -d" " -f2,3 tmp`basename $0`3.txt > tmp`basename $0`4.txt

# plot
PlotIbd.R tmp`basename $0`4.txt $3
