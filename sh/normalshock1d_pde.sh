#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/normalshock1d_pde

rm ../lp_output/normalshock1d_pde.txt

./lp -i input/input_normalshock1d_pde -o ../lp_output/normalshock1d_pde > ../lp_output/normalshock1d_pde.txt
