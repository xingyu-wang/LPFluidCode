#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_rayleightaylor_SPH_coarse

rm ../lp_output/review/rayleightaylor_SPH_coarse.txt

./lp -i input/review/input_rayleightaylor_lw_nosplit_SPH_coarse -o ../lp_output/review/lp_rayleightaylor_SPH_coarse > ../lp_output/review/rayleightaylor_SPH_coarse.txt
