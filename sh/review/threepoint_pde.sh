#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_threepoint_pde

rm ../lp_output/review/threepoint_pde.txt

./lp -i input/review/input_threepoint2d_pde -o ../lp_output/review/lp_threepoint_pde > ../lp_output/review/threepoint_pde.txt
