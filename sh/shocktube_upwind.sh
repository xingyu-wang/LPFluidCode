#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_shocktube_upwind_01

rm ../lp_output/shocktube_upwind_01.txt

./lp -i input/input_shocktube2d_solidb_upwind -o ../lp_output/lp_shocktube_upwind_01 > ../lp_output/shocktube_upwind_01.txt
