#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/performance/lp_kh

rm ../lp_output/performance/kh.txt

./lp -i input/input_kelvinhelmholtz2d_pde_24 -o ../lp_output/performance/lp_kh > ../lp_output/performance/kh.txt
