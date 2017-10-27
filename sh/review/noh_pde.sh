#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_noh_pde

rm ../lp_output/review/noh_pde.txt

./lp -i input/review/input_noh2d_pde -o ../lp_output/review/lp_noh_pde > ../lp_output/review/noh_pde.txt
