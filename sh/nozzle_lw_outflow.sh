#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_lw_outflow

rm ../lp_output/nozzle_lw_outflow.txt

./lp -i input/input_nozzle_lw_outflow -o ../lp_output/lp_nozzle_lw_outflow > ../lp_output/nozzle_lw_outflow.txt
