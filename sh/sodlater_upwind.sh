#!/bin/bash
#PBS -l nodes=1:ppn=20

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_sodshocktubelater_upwind

rm ../lp_output/sodshocktubelater_upwind.txt

./lp -i input/input_sodshocktubelater_upwind -o ../lp_output/lp_sodshocktubelater_upwind > ../lp_output/sodshocktubelater_upwind.txt
