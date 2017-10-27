#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_dambreak_lw_nosplit_pde_100_r

rm ../lp_output/dambreak_lw_nosplit_pde_100_r.txt

./lp -i input/input_dambreak_lw_nosplit_100 -o ../lp_output/lp_dambreak_lw_nosplit_pde_100_r > ../lp_output/dambreak_lw_nosplit_pde_100_r.txt
