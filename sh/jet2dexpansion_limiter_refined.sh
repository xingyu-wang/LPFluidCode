#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dexpansion_limiter_refined

rm ../lp_output/jet2dexpansion_limiter_refined.txt

./lp -i input/input_jet2xpansion_limiter_refined -o ../lp_output/jet2dexpansion_limiter_refined> ../lp_output/jet2dexpansion_limiter_refined.txt
