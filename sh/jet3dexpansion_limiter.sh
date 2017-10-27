#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet3dexpansion_limiter

rm ../lp_output/jet3dexpansion_limiter.txt

./lp -i input/input_jet3dexpansion_limiter -o ../lp_output/jet3dexpansion_limiter > ../lp_output/jet3dexpansion_limiter.txt
