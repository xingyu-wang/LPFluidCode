#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktube3d_upwind_test

rm ../lp_output/sodshocktube3d_upwind_test.txt

./lp -i input/input_sodshocktube3d_upwind -o ../lp_output/lp_sodshocktube3d_upwind_test > ../lp_output/sodshocktube3d_upwind_test.txt
