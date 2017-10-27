#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dmerge_coarse_quadrant3

rm ../lp_output/jet2dmerge_coarse_quadrant3.txt

./lp -i input/input_jet2dmerge_coarse -o ../lp_output/jet2dmerge_coarse_quadrant3 > ../lp_output/jet2dmerge_coarse_quadrant3.txt
