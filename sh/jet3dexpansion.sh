#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet3dexpansion

rm ../lp_output/jet3dexpansion.txt

./lp -i input/input_jet3dexpansion -o ../lp_output/jet3dexpansion > ../lp_output/jet3dexpansion.txt
