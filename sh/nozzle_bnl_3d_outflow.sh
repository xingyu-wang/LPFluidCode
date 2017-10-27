#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_bnl_3d_outflow

rm ../lp_output/nozzle_bnl_3d_outflow.txt

./lp -i input/input_nozzle_bnl_3d_outflow -o ../lp_output/lp_nozzle_bnl_3d_outflow > ../lp_output/nozzle_bnl_3d_outflow.txt
