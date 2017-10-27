#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_stiff_free_surface_lw

rm ../lp_output/stiff_free_surface_lw.txt

./lp -i input/input_stiff_free_surface_lw -o ../lp_output/lp_stiff_free_surface_lw > ../lp_output/stiff_free_surface_lw.txt
