#!/bin/bash
#SBATCH --job-name=clusterKeepAlive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.Out
#SBATCH --partition=hpc

set -euo pipefail

echo "[$(date '+%Y-%m-%d %H:%M:%S')] keep-alive job started on ${SLURM_JOB_NODELIST:-unknown}"
echo "This job intentionally performs minimal work to keep one allocation active."

# Keep one process alive with negligible load.
while true; do
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] heartbeat"
  sleep 300
done
