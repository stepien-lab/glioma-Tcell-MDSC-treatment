#!/bin/sh
#SBATCH --job-name=OCHiPerGator       # Use the same name as the MATLAB program\
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)\
#SBATCH --mail-user=hannahanderson@ufl.edu    # Where to send mail\
#SBATCH --account=tstepien
#SBATCH --qos=tstepien-b
#SBATCH --nodes=1                       # Use one node\
#SBATCH --ntasks=1                      # Run a single task\
#SBATCH --cpus-per-task=4              # Check QoS/queue to see how many CPUS are available\
#SBATCH --mem-per-cpu=2gb               # Memory per processor - max 2gb for MATLAB\
#SBATCH --time=00:30:00                 # Time limit hrs:min:sec\ % takes about 30 hours for 10,000 parameter sets
#SBATCH --output=logs/MATLAB_%j.txt     # Output and error log - create "logs" folder first before running SLURM script\
pwd; hostname; date

module load matlab/2019a
./OCHiPerGator

date
exit
