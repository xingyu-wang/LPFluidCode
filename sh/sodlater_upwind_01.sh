#!/bin/bash
#PBS -l nodes=1:ppn=20

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_sodshocktubelater_upwind_01

#rm ../lp_output/sodshocktubelater_upwind_01.txt

./lp -i input/input_sodshocktubelater_upwind_01 -o ../lp_output/lp_sodshocktubelater_upwind_01 > ../lp_output/sodshocktubelater_upwind_01.txt
