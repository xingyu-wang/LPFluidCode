#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_bnl_test

rm ../lp_output/nozzle_bnl_test.txt

./lp -i input/input_nozzle_bnl_3d_outflow -o ../lp_output/lp_nozzle_bnl_test > ../lp_output/nozzle_bnl_test.txt
