#!/bin/bash
#PBS -l nodes=vn04:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktube_upwind_old

rm ../lp_output/sodshocktube_upwind_old.txt

./lp -i input/input_sodshocktube_upwind -o ../lp_output/lp_sodshocktube_upwind_old > ../lp_output/sodshocktube_upwind_old.txt
