#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_nozzle

rm ../lp_output/nozzle.txt

./lp -i input/input_nozzle -o ../lp_output/lp_nozzle > ../lp_output/nozzle.txt
