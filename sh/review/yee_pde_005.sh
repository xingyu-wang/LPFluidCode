#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_yee_pde_random_005

rm ../lp_output/review/yee_pde_random_005.txt

./lp -i input/review/input_yee2d_pde_005 -o ../lp_output/review/lp_yee_pde_random_005 > ../lp_output/review/yee_pde_random_005.txt
