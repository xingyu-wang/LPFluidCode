#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_shocktube_lw_nosplit

rm ../lp_output/shocktube_lw_nosplit.txt

./lp -i input/input_shocktube2d_solidb_lw_nosplit -o ../lp_output/lp_shocktube_lw_nosplit > ../lp_output/shocktube_lw_nosplit.txt
