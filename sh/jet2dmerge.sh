#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dmerge

rm ../lp_output/jet2dmerge.txt

./lp -i input/input_jet2dmerge -o ../lp_output/jet2dmerge > ../lp_output/jet2dmerge.txt
