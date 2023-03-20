#! usr/bin/bash

# i=1
# name=Tux
# all=
# echo $i $name

for i in {1..20}
do
	inst=G$i
	echo $inst
    ./an_ss_ge_fi_vdeg_omp -v -l example/$inst.txt -r 500 -s 10000 -b0 0.1 -b1 3 -sched exp > result/$inst.txt
done
