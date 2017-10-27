#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_gresho_upwind

rm ../lp_output/gresho_upwind.txt

./lp -i input/input_gresho2d_solidb -o ../lp_output/lp_gresho_upwind > ../lp_output/gresho_upwind.txt
