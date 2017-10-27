#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/simplewave1d_upwind

rm ../lp_output/simplewave1d_upwind.txt

./lp -i input/input_simplewave1d_solidb_upwind -o ../lp_output/simplewave1d_upwind > ../lp_output/simplewave1d_upwind.txt
