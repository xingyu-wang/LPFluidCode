#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_ballrotate_upwind_sph

rm ../lp_output/ballrotate_upwind_sph.txt

./lp -i input/input_ballrotate_upwind_sph -o ../lp_output/lp_ballrotate_upwind_sph > ../lp_output/ballrotate_upwind_sph.txt
