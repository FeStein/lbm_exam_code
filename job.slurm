#!/bin/bash
#SBATCH --mail-type=End
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err            # to report error
#SBATCH -p mv-ltm               # ltm partition
#SBATCH --ntasks=16              # Number of processors
#SBATCH -t 3000                  #runtime in min
#SBATCH --mem-per-cpu=2GB

for U in 0.5 1 1.5
do
    FOLDERNAME="${U}_U"
    mkdir ${FOLDERNAME}
    ./channelFlow ${U}
    mv *.vtk *.log ${FOLDERNAME} 
done


