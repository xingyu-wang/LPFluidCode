#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dexpansion

rm ../lp_output/jet2dexpansion.txt

./lp -i input/input_jet2dexpansion -o ../lp_output/jet2dexpansion > ../lp_output/jet2dexpansion.txt
