#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_rayleightaylor_SPH_fine

rm ../lp_output/review/rayleightaylor_SPH_fine.txt

./lp -i input/review/input_rayleightaylor_lw_nosplit_SPH_fine -o ../lp_output/review/lp_rayleightaylor_SPH_fine > ../lp_output/review/rayleightaylor_SPH_fine.txt
