#!/bin/bash
#SBATCH -J unscheduled_gps_paper_pipeline
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake,cclake-himem,icelake,icelake-himem,skylake,skylake-himem
#SBATCH -o unscheduled_gps_paper_pipeline/%j.out

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load r-4.0.2-gcc-5.4.0-xyx46xb
export R_LIBS=/home/tw395/R/4.0.2/libs:$R_LIBS

workdir="$SLURM_SUBMIT_DIR"

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on mcadr node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > logs/machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat logs/machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n\n"

# Needed to run the code written to .bashrc by the conda init program
source  /home/tw395/.bashrc
# Needed to use my libs
source  /home/tw395/.bash_profile

conda activate gps_paper_pipeline

snakemake --nolock --cores 6 -j 300 --keep-going --use-conda --scheduler greedy --rerun-incomplete "${@}"
