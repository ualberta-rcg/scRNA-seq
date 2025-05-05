#!/bin/bash
#SBATCH --job-name=rstudio
#SBATCH --account=def-xxxxxx		# Adjust the account
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         # Adjust CPU resources as needed
#SBATCH --mem=8G                 # Adjust memory as needed
#SBATCH --time=04:00:00           # Set the job runtime limit
#SBATCH --output=job_%j.log   # Log file
#SBATCH -e job_%j.err

# Load Apptainer module (if needed)
module load apptainer

# Define ports
RPORT=$((RANDOM % 10000 + 10000))       # Default RStudio Server port
SSHPORT=$((RANDOM % 20000 + 20000))     # Random number from 20000-40000

# Work directory for Apptainer
WORKDIR="/$HOME/rstudio/job_$SLURM_JOB_ID"
mkdir -p $WORKDIR
mkdir $WORKDIR/lib
mkdir $WORKDIR/run
mkdir $WORKDIR/home

# Get compute node hostname
NODE_HOSTNAME=$(hostname)

echo "Job running on: $NODE_HOSTNAME"
echo "Setting up SSH tunnel for RStudio..."

# Show access instructions
echo -e "\nTo connect to RStudio from your local machine, run:"
echo -e "ssh -L $SSHPORT:${NODE_HOSTNAME}:$RPORT $USER@cedar.alliancecan.ca\n"
echo "Then, open your web browser and go to:"
echo -e "http://localhost:$SSHPORT\n"


# Start rstudio
apptainer exec --workdir $WORKDIR --home $WORKDIR/home --bind $WORKDIR/lib:/var/lib/rstudio-server --bind $WORKDIR/run:/var/run/rstudio-server /$HOME/scRNA-seq.sif rserver --www-port=$RPORT --server-daemonize=0 --server-user=$(whoami)

