#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_gresho_pde_random_plotneighbour

rm ../lp_output/review/gresho_pde_random_plotneighbour.txt

./lp -i input/review/input_gresho2d_pde -o ../lp_output/review/lp_gresho_pde_random_plotneighbour > ../lp_output/review/gresho_pde_random_plotneighbour.txt
