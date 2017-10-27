#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2d

rm ../lp_output/jet2d.txt

./lp -i input/input_jet2d -o ../lp_output/jet2d > ../lp_output/jet2d.txt
