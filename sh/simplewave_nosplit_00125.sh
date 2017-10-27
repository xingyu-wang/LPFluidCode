#!/bin/bash
#PBS -l nodes=1:ppn=20

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_simplewave_nosplit_00125

rm ../lp_output/simplewave_nosplit_00125.txt

./lp -i input/input_simplewave2d_solidb_lw_nosplit_00125 -o ../lp_output/lp_simplewave_nosplit_00125 > ../lp_output/simplewave_nosplit_00125.txt
