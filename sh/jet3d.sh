#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet3d

rm ../lp_output/jet3d.txt

./lp -i input/input_jet3d -o ../lp_output/jet3d > ../lp_output/jet3d.txt
