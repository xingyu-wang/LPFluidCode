#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_simplewave_nosplit

rm ../lp_output/simplewave_nosplit.txt

./lp -i input/input_simplewave2d_solidb_lw_nosplit -o ../lp_output/lp_simplewave_nosplit > ../lp_output/simplewave_nosplit.txt
