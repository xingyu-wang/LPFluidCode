#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_gresho_lw_nosplit_cfl02

#rm ../lp_output/gresho_lw_nosplit_cfl02.txt

./lp -i input/input_gresho2d_solidb_lw_nosplit_cfl02 -o ../lp_output/lp_gresho_lw_nosplit_cfl02 > ../lp_output/gresho_lw_nosplit_cfl02.txt
