#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_kh_pde_random_new

rm ../lp_output/review/kh_pde_random_new.txt

./lp -i input/review/input_kelvinhelmholtz2d_pde -o ../lp_output/review/lp_kh_pde_random_new > ../lp_output/review/kh_pde_random_new.txt
