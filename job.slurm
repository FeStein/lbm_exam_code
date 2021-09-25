#!/bin/bash
#SBATCH --mail-type=End
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err            # to report error
#SBATCH -p mv-ltm               # ltm partition
#SBATCH --ntasks=16              # Number of processors
#SBATCH -t 500                  #runtime in min
#SBATCH --mem-per-cpu=2GB

for U in 1 2 5 10 
do
    FOLDERNAME="${U}_dt"
    mkdir ${FOLDERNAME}
    ./channelFlow ${U}
    mv *.vtk *.log ${FOLDERNAME} 
done


