#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/performance/lp_kh_2

#rm ../lp_output/performance/kh_2.txt

./lp -i input/input_kelvinhelmholtz2d_pde_2 -o ../lp_output/performance/lp_kh_2 > ../lp_output/performance/kh_2.txt
