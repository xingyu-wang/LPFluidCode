#PBS -N interpolation
#PBS -l nodes=1:ppn=24
#PBS -l walltime=192:00:00

cd /home/xingyu/src_08_12_15/data_processor

./interpolation > screen1.txt
