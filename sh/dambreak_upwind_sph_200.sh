#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_dambreak_upwind_sph_200

rm ../lp_output/dambreak_upwind_sph_200.txt

./lp -i input/input_dambreak_upwind_200 -o ../lp_output/lp_dambreak_upwind_sph_200 > ../lp_output/dambreak_upwind_sph_200.txt
