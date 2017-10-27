#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/normalshock1d_sph

rm ../lp_output/normalshock1d_sph.txt

./lp -i input/input_normalshock1d_sph -o ../lp_output/normalshock1d_sph > ../lp_output/normalshock1d_sph.txt
