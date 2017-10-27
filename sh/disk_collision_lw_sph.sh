#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_disk_collision_lw_sph

rm ../lp_output/disk_collision_lw_sph.txt

./lp -i input/input_disk_collision_lw_sph -o ../lp_output/lp_disk_collision_lw_sph > ../lp_output/disk_collision_lw_sph.txt
