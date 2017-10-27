#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_dambreak_lw_nosplit_sph

rm ../lp_output/dambreak_lw_nosplit_sph.txt

./lp -i input/input_dambreak_lw_nosplit_sph -o ../lp_output/lp_dambreak_lw_nosplit_sph > ../lp_output/dambreak_lw_nosplit_sph.txt
