#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_threepoint_sph_v

rm ../lp_output/review/threepoint_sph_v.txt

./lp -i input/review/input_threepoint2d_sph -o ../lp_output/review/lp_threepoint_sph_v > ../lp_output/review/threepoint_sph_v.txt
