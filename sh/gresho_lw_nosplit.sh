#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_gresho_lw_nosplit

rm ../lp_output/gresho_lw_nosplit.txt

./lp -i input/input_gresho2d_solidb_lw_nosplit -o ../lp_output/lp_gresho_lw_nosplit > ../lp_output/gresho_lw_nosplit.txt
