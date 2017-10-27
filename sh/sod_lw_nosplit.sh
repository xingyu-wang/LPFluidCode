#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktube_lw_nosplit

rm ../lp_output/sodshocktube_lw_nosplit.txt

./lp -i input/input_sodshocktube_lw_nosplit -o ../lp_output/lp_sodshocktube_lw_nosplit > ../lp_output/sodshocktube_lw_nosplit.txt
