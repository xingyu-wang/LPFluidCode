#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_gresho_test

rm ../lp_output/review/gresho_test.txt

./lp -i input/review/input_gresho2d_pde -o ../lp_output/review/lp_gresho_test > ../lp_output/review/gresho_test.txt
