#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_rothe_delete

rm ../lp_output/nozzle_rothe_delete.txt

./lp -i input/input_nozzle_rothe_delete -o ../lp_output/lp_nozzle_rothe_delete > ../lp_output/nozzle_rothe_delete.txt
