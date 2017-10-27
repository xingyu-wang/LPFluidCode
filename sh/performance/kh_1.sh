#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/performance/lp_kh_1

#rm ../lp_output/performance/kh_1.txt

./lp -i input/input_kelvinhelmholtz2d_pde_1 -o ../lp_output/performance/lp_kh_1 > ../lp_output/performance/kh_1.txt
