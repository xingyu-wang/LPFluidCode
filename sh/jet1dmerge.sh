#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/jet1dmerge

rm ../lp_output/jet1dmerge.txt

./lp -i input/input_jet1dmerge -o ../lp_output/jet1dmerge > ../lp_output/jet1dmerge.txt
