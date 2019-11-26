#!/usr/bin/env bash
#PBS -q expressbw
#PBS -P x52
#PBS -l walltime=04:00:00
#PBS -l ncpus=28
#PBS -l mem=44044mb
#PBS -N ML_5_N
##PBS -wd
#PBS -joe

#PBS -v count,max
##PBS -v 1,5

# Test bash script to demonstrate resubmissions
script_name='resub.sh'

#script_name='./diablo > output.dat &'


# Set default values of count and max
if [ -z $count ]; then
    count=1
fi

if [ -z $max ]; then
    max=$count
fi

# Log submission counters
echo "Run $count of $max"

if [ $count -gt 1 ]; then
    # Copy previous restart files to input path
    cd $PBS_O_WORKDIR
    ./copy_start_saved 
    echo 'Copy files here'
fi

# Run the model
mpirun -n 16 ./diablo > output.dat

((count++))

if [ $count -le $max ]; then
    echo "Resubmitting model"
    cd $PBS_O_WORKDIR
    qsub -v count=$count,max=$max $script_name
else
    echo "Last submission"
fi
