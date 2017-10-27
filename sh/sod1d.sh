#!/bin/bash
#PBS -l nodes=1:ppn=4

cd /home/xingyu/Downloads/lp_code/src_08_12_15

#rm -r ../lp_output/sod1d_new_sph

#rm ../lp_output/sod1d_new_sph.txt

./lp -i input/input_sodshocktube1d -o ../lp_output/sod1d_new_pde > ../lp_output/sod1d_new_pde.txt
