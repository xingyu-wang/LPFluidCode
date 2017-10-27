#!/bin/bash
#PBS -l nodes=1:ppn=20

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktube_upwind_pde_0025

rm ../lp_output/sodshocktube_upwind_pde_0025.txt

./lp -i input/input_sodshocktube_upwind_pde_0025 -o ../lp_output/lp_sodshocktube_upwind_pde_0025 > ../lp_output/sodshocktube_upwind_pde_0025.txt
