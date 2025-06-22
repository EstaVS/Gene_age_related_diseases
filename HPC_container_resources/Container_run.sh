#!/bin/bash

#SBATCH --job-name=ops-rstudio
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=16G            # Increased from 2G → 16G
#SBATCH --cpus-per-task=4     # Increased from 1 → 4 CPUs
#SBATCH --signal=USR2
#SBATCH --time=24:00:0

# Set your container path
CONTAINER="/scratch/prj/bmb_phyl_chron_disease/esta_proj/bioconductor_latest.sif"

# Verify container exists
if [ ! -f "$CONTAINER" ]; then
    echo "Error: Container not found at $CONTAINER" >&2
    exit 1
fi

# Get connection info
export PASSWORD=$(openssl rand -base64 15)
PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')

cat 1>&2 <<END
1. SSH tunnel from your workstation:
   ssh -i ~/.ssh/create_msc -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

2. Connect to: http://localhost:8787
   Username: ${USER}
   Password: ${PASSWORD}
END

# Create temp directory in your scratch space
LOCAL_TMP="/scratch/users/k24058218/rstudio_temp_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_TMP}"/{var-rstudio-server,data-rstudio-server} || {
    echo "Error: Failed to create temp directory at ${LOCAL_TMP}" >&2
    exit 1
}

# Launch RStudio in container
singularity exec \
  -B "${LOCAL_TMP}:/tmp/rstudio-tmp" \
  -B "${HOME}" \
  -B "/scratch/users/k24058218" \
  -B "/scratch/prj/bmb_phyl_chron_disease/esta_proj:/scratch/prj/bmb_phyl_chron_disease/esta_proj" \
  --env "TMPDIR=/tmp/rstudio-tmp" \
  --env "PASSWORD=${PASSWORD}" \
  "$CONTAINER" \
  bash -c "
    echo 'directory=/tmp/rstudio-tmp/var-rstudio-server' > /tmp/rstudio-tmp/database.conf
    rserver \
      --server-user ${USER} \
      --www-port ${PORT} \
      --server-data-dir /tmp/rstudio-tmp/data-rstudio-server \
      --secure-cookie-key-file /tmp/rstudio-tmp/data-rstudio-server/secure-cookie-key \
      --database-config-file=/tmp/rstudio-tmp/database.conf \
      --auth-none=0 \
      --auth-pam-helper-path=pam-helper
  " || {
    echo "RStudio Server failed with status $?" >&2
  }

# Cleanup (only if successful)
[ $? -eq 0 ] && rm -rf "${LOCAL_TMP}"
