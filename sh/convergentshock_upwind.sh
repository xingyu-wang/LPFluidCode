#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_convergentshock_upwind

rm ../lp_output/convergentshock_upwind.txt

./lp -i input/input_convergentshock2d_solidb_upwind -o ../lp_output/lp_convergentshock_upwind > ../lp_output/convergentshock_upwind.txt
