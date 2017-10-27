#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/review/lp_sedov_sph

rm ../lp_output/review/sedov_sph.txt

./lp -i input/review/input_sedov2d_sph -o ../lp_output/review/lp_sedov_sph > ../lp_output/review/sedov_sph.txt
