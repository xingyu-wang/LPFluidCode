#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/review/lp_yee_sph_norandom_04

#rm ../lp_output/review/yee_sph_norandom_04.txt

./lp -i input/review/input_yee2d_sph_04 -o ../lp_output/review/lp_yee_sph_norandom_04 > ../lp_output/review/yee_sph_norandom_04.txt
