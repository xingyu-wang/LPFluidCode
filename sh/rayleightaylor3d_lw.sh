#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_rayleightaylor3d_lw

#rm ../lp_output/rayleightaylor3d_lw.txt

./lp -i input/input_rayleightaylor3d_lw_nosplit -o ../lp_output/lp_rayleightaylor3d_lw_random > ../lp_output/rayleightaylor3d_lw_random.txt
