#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

rm -r ../lp_output/jet1d

rm ../lp_output/jet1d.txt

./lp -i input/input_jet1d -o ../lp_output/jet1d > ../lp_output/jet1d.txt
