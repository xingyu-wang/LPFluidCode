#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dexpansion_limiter

rm ../lp_output/jet2dexpansion_limiter.txt

./lp -i input/input_jet2dexpansion_limiter -o ../lp_output/jet2dexpansion_limiter> ../lp_output/jet2dexpansion_limiter.txt
