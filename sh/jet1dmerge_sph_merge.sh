#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/jet1dmerge_sph_merge

rm ../lp_output/jet1dmerge_sph_merge.txt

./lp -i input/input_jet1dmerge_sph_merge -o ../lp_output/jet1dmerge_sph_merge > ../lp_output/jet1dmerge_sph_merge.txt
