{\rtf1\ansi\ansicpg1252\cocoartf2708
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/sh\
#SBATCH --job-name=OCHiPerGator       # Use the same name as the MATLAB program\
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)\
#SBATCH --mail-user=hannahanderson@ufl.edu    # Where to send mail\
#SBATCH --account=hannahanderson\
#SBATCH --qos=tstepien-b\
#SBATCH --nodes=1                       # Use one node\
#SBATCH --ntasks=1                      # Run a single task\
#SBATCH --cpus-per-task=2              # Check QoS/queue to see how many CPUS are available\
#SBATCH --mem-per-cpu=2gb               # Memory per processor - max 2gb for MATLAB\
#SBATCH --time=00:30:00                 # Time limit hrs:min:sec\
#SBATCH --output=logs/MATLAB_%j.txt     # Output and error log - create "logs" folder first before running SLURM script\
pwd; hostname; date\
\
module load matlab\
./OCHiPerGator\
\
date\
exit}