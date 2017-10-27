#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_rayleightaylor_upwind

#rm ../lp_output/rayleightaylor_upwind.txt

./lp -i input/input_rayleightaylor_upwind -o ../lp_output/lp_rayleightaylor_upwind_refine > ../lp_output/rayleightaylor_upwind_refine.txt
