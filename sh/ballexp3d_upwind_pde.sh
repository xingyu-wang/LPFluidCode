#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_ballexp3d_upwind_pde

rm ../lp_output/ballexp3d_upwind_pde.txt

./lp -i input/input_ballexp3d_upwind_pde -o ../lp_output/lp_ballexp3d_upwind_pde > ../lp_output/ballexp3d_upwind_pde.txt
