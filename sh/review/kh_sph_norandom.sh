#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_kh_sph_norandom

rm ../lp_output/review/kh_sph_norandom.txt

./lp -i input/review/input_kelvinhelmholtz2d_sph -o ../lp_output/review/lp_kh_sph_norandom > ../lp_output/review/kh_sph_norandom.txt
