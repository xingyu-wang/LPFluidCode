#PBS -N ballexp
#PBS -l nodes=1:towel:ppn=12 
#PBS -V 
cd $PBS_O_WORKDIR 
make clean
rm -r out
make -f makefile_vogon
./lp -i ./input/input_ballexp3d -o out > log

