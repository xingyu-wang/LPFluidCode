#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_rayleightaylor_lw_label

#rm ../lp_output/rayleightaylor_lw_voronoi_test.txt

./lp -i input/input_rayleightaylor_lw_nosplit -o ../lp_output/lp_rayleightaylor_lw_refine_label_random1 > ../lp_output/rayleightaylor_lw_refine_label_random1.txt
