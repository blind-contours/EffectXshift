#!/bin/bash
# Job name:
#SBATCH --job-name=4_c_quad_mediation_EffectXshift
#
# Partition:
#SBATCH --partition=savio3
#

#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat

# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --exclusive
#
# Wall clock limit:
#SBATCH --time 72:00:00
#
## Command(s) to run (example):
module load r/4.0.3

### Run Simulation
R CMD BATCH --no-save ../03_run_med_simulation_cont_quad_4.R EffectXshift_mediation_cont_quad_4.Rout
