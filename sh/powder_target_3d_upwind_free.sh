#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

rm -r ../lp_output/lp_powder_target_3d_upwind_free

rm ../lp_output/powder_target_3d_upwind_free.txt

./lp -i input/input_powder_target_3d_upwind_free -o ../lp_output/lp_powder_target_3d_upwind_free > ../lp_output/powder_target_3d_upwind_free.txt
