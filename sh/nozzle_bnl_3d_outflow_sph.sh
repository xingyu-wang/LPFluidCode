#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/src_08_12_15

#rm -r ../lp_output/lp_nozzle_bnl_3d_outflow_sph

#rm ../lp_output/nozzle_bnl_3d_outflow_sph.txt

./lp -i input/input_nozzle_bnl_3d_outflow_sph -o ../lp_output/lp_nozzle_bnl_3d_outflow_sph_long > ../lp_output/nozzle_bnl_3d_outflow_sph_long.txt
