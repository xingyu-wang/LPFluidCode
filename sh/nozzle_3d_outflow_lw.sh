#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_3d_outflow_lw

rm ../lp_output/nozzle_3d_outflow_lw.txt

./lp -i input/input_nozzle_3d_outflow_lw -o ../lp_output/lp_nozzle_3d_outflow_lw > ../lp_output/nozzle_3d_outflow_lw.txt
