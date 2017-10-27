#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_ballrotate_lw_pde

rm ../lp_output/ballrotate_lw_pde.txt

./lp -i input/input_ballrotate_lw_pde -o ../lp_output/lp_ballrotate_lw_pde > ../lp_output/ballrotate_lw_pde.txt
