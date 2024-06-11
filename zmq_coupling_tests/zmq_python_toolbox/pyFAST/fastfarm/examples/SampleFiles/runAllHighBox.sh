#!/bin/bash
#SBATCH --job-name=highBox
#SBATCH --output log.highBox
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=36
#SBATCH --time=4:00:00
#SBATCH --mem=150G
#SBATCH --account=osw

source $HOME/.bash_profile

echo "Working directory is" $SLURM_SUBMIT_DIR
echo "Job name is" $SLURM_JOB_NAME
echo "Job ID is " $SLURM_JOBID
echo "Job took the following nodes (SLURM_NODELIST)" $SLURM_NODELIST
echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "Starting job at: " $(date)

module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles
module load mkl/2020.1.217
module load comp-intel/2020.1.217

# ********************************** USER INPUT ********************************** #
turbsimbin='/home/rthedin/local/local_openfast_intelCompilers/bin/turbsim'
basepath='/projects/shellwind/rthedin/Task2_2regis'

condList=('Cond00_v08.6_PL0.2_TI10' 'Cond01_v10.6_PL0.2_TI10' 'Cond02_v12.6_PL0.2_TI10')

caseList=('Case00_wdirp00_WSfalse_YMfalse' 'Case01_wdirp00_WStrue_YMfalse')

nSeeds=6
nTurbines=12
# ******************************************************************************** #

rampercpu=$((149000/36))

for cond in ${condList[@]}; do
    for case in ${caseList[@]}; do
        for ((seed=0; seed<$nSeeds; seed++)); do
            for ((t=1; t<=$nTurbines; t++)); do
                dir=$(printf "%s/%s/%s/Seed_%01d/TurbSim" $basepath $cond $case $seed)
                echo "Submitting $dir/HighT$t.inp"
                srun -n1 -N1 --exclusive --mem-per-cpu=$rampercpu $turbsimbin $dir/HighT$t.inp > $dir/log.hight$t.seed$seed.txt 2>&1 &
                sleep 0.1
            done
        done
    done
done

wait

echo "Ending job at: " $(date)
