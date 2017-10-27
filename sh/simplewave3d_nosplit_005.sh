#!/bin/bash
#PBS -l nodes=1:ppn=20

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_simplewave3d_nosplit_005

#rm ../lp_output/simplewave3d_nosplit_005.txt

./lp -i input/input_simplewave3d_solidb_lw_nosplit_005 -o ../lp_output/lp_simplewave3d_nosplit_sph_005 > ../lp_output/simplewave3d_nosplit_sph_005.txt
