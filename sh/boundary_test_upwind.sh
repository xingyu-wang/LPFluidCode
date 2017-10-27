#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_boundary_test_upwind

rm ../lp_output/boundary_test_upwind.txt

./lp -i input/input_boundary_test_upwind -o ../lp_output/lp_boundary_test_upwind > ../lp_output/boundary_test_upwind.txt
