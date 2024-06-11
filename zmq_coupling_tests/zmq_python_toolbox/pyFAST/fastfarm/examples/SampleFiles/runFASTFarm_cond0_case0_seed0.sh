#!/bin/bash
#SBATCH --job-name=runFF
#SBATCH --output log.fastfarm_c0_c0_seed0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=10-00
#SBATCH --account=mmc

source $HOME/.bash_profile

echo "Working directory is" $SLURM_SUBMIT_DIR
echo "Job name is" $SLURM_JOB_NAME
echo "Job ID is " $SLURM_JOBID
echo "Job took the following nodes (SLURM_NODELIST)" $SLURM_NODELIST
echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "Starting job at: " $(date)

nodelist=`scontrol show hostname $SLURM_NODELIST`
nodelist=($nodelist)
echo "Formatted list of nodes is: $nodelist"

module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles
module load mkl/2020.1.217
module load comp-intel/2020.1.217

# ********************************** USER INPUT ********************************** #
fastfarmbin='/home/rthedin/local/local_openfast_intelCompilers_3.1.0_openmpPrints/bin/FAST.Farm'
basepath='/projects/shellwind/rthedin/Task2_2regis'

cond='Cond00_v08.6_PL0.2_TI10'
case='Case00_wdirp00_WSfalse_YMfalse_12fED_12ADyn'

seed=0
# ******************************************************************************** #

dir=$(printf "%s/%s/%s/Seed_%01d" $basepath $cond $case $seed)
cd $dir
echo "Submitting $dir/FFarm_mod.fstf"
$fastfarmbin $dir/FFarm_mod.fstf > $dir/log.fastfarm.seed$seed.txt 2>&1

echo "Ending job at: " $(date)
