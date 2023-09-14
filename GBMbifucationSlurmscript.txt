#!/bin/sh
#SBATCH --job-name=GBMbifurcationanalysis        # Use the same name as the MATLAB program
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hannahanderson@ufl.edu    # Where to send mail
#SBATCH --account=hannahanderson
#SBATCH --qos=tstepien-b
#SBATCH --nodes=1                       # Use one node
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=1               # Check QoS/queue to see how many CPUS are available
#SBATCH --mem-per-cpu=2gb               # Memory per processor - max 2gb for MATLAB
#SBATCH --time=00:10:00                 # Time limit hrs:min:sec
#SBATCH --output=logs/MATLAB_%j.txt     # Output and error log - create "logs" folder first before running SLURM script
pwd; hostname; date

module load matlab
./GBMbifurcationanalysis

date
exit