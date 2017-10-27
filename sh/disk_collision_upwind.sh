#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_disk_collision_upwind

rm ../lp_output/disk_collision_upwind.txt

./lp -i input/input_disk_collision_upwind -o ../lp_output/lp_disk_collision_upwind > ../lp_output/disk_collision_upwind.txt
