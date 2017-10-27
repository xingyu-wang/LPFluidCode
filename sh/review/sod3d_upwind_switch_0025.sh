#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_sodshocktube3d_upwind_switch_0025

rm ../lp_output/review/sodshocktube3d_upwind_switch_0025.txt

./lp -i input/review/input_sodshocktube3d_upwind_switch_0025 -o ../lp_output/review/lp_sodshocktube3d_upwind_switch_0025 > ../lp_output/review/sodshocktube3d_upwind_switch_0025.txt
