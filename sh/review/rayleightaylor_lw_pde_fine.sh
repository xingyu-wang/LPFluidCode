#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_rayleightaylor_PDE_fine

rm ../lp_output/review/rayleightaylor_PDE_fine.txt

./lp -i input/review/input_rayleightaylor_lw_nosplit_PDE_fine -o ../lp_output/review/lp_rayleightaylor_PDE_fine > ../lp_output/review/rayleightaylor_PDE_fine.txt
