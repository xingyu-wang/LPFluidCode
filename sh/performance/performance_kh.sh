#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/performance/lp_kh_*;
rm ../lp_output/performance/kh_*.txt;

./lp -i input/input_kelvinhelmholtz2d_pde_24 -o ../lp_output/performance/lp_kh_24 > ../lp_output/performance/kh_24.txt

./lp -i input/input_kelvinhelmholtz2d_pde_16 -o ../lp_output/performance/lp_kh_16 > ../lp_output/performance/kh_16.txt

./lp -i input/input_kelvinhelmholtz2d_pde_8 -o ../lp_output/performance/lp_kh_8 > ../lp_output/performance/kh_8.txt

./lp -i input/input_kelvinhelmholtz2d_pde_4 -o ../lp_output/performance/lp_kh_4 > ../lp_output/performance/kh_4.txt

./lp -i input/input_kelvinhelmholtz2d_pde_2 -o ../lp_output/performance/lp_kh_2 > ../lp_output/performance/kh_2.txt

./lp -i input/input_kelvinhelmholtz2d_pde_1 -o ../lp_output/performance/lp_kh_1 > ../lp_output/performance/kh_1.txt
