#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_sodshocktube_upwind_switch_005

rm ../lp_output/review/sodshocktube_upwind_switch_005.txt

./lp -i input/review/input_sodshocktube_upwind_switch_005 -o ../lp_output/review/lp_sodshocktube_upwind_switch_005 > ../lp_output/review/sodshocktube_upwind_switch_005.txt
