#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_bnl

rm ../lp_output/nozzle_bnl.txt

./lp -i input/input_nozzle_bnl -o ../lp_output/lp_nozzle_bnl > ../lp_output/nozzle_bnl.txt
