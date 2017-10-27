#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dmerge_first_soarse

rm ../lp_output/jet2dmerge_first_coarse.txt

./lp -i input/input_jet2dmerge_first_coarse -o ../lp_output/jet2dmerge_first_coarse > ../lp_output/jet2dmerge_first_coarse.txt
