#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/jet2dmerge_first

rm ../lp_output/jet2dmerge_first.txt

./lp -i input/input_jet2dmerge_first -o ../lp_output/jet2dmerge_first > ../lp_output/jet2dmerge_first.txt
