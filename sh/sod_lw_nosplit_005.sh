#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktube_lw_nosplit_005

rm ../lp_output/sodshocktube_lw_nosplit_005.txt

./lp -i input/input_sodshocktube_lw_nosplit_005 -o ../lp_output/lp_sodshocktube_lw_nosplit_005 > ../lp_output/sodshocktube_lw_nosplit_005.txt
