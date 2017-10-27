#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/performance/lp_kh_8

#rm ../lp_output/performance/kh_8.txt

./lp -i input/input_kelvinhelmholtz2d_pde_8 -o ../lp_output/performance/lp_kh_8 > ../lp_output/performance/kh_8.txt
