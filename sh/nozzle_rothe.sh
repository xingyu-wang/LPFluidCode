#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle_rothe

rm ../lp_output/nozzle_rothe.txt

./lp -i input/input_nozzle_rothe -o ../lp_output/lp_nozzle_rothe > ../lp_output/nozzle_rothe.txt
