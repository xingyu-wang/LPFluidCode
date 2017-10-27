#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_kh

rm ../lp_output/kh.txt

./lp -i input/input_kelvinhelmholtz2d_pde -o ../lp_output/lp_kh > ../lp_output/kh.txt
