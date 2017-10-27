#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_dambreak_upwind

rm ../lp_output/dambreak_upwind.txt

./lp -i input/input_dambreak_upwind -o ../lp_output/lp_dambreak_upwind > ../lp_output/dambreak_upwind.txt
