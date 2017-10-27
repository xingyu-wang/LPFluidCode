#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_kh_switch_compare_random

rm ../lp_output/review/kh_switch_compare_random.txt

./lp -i input/review/input_kelvinhelmholtz2d_switch -o ../lp_output/review/lp_kh_switch_compare_random > ../lp_output/review/kh_switch_compare_random.txt
