#!/bin/bash

#$ -N SuperFinal # Job name
#$ -t 1-1500     # Number of jobs
#$ -q all.q      # Queue. Use long.q for run time >8h 
#$ -tc 128       # Max concurrent jobs
#$ -cwd          # Run in current directory
#$ -o output/    # Direct output to subdirectory
#$ -e output/    # Direct output to subdirectory

R CMD BATCH job.R output/$JOB_NAME-$SGE_TASK_ID.Rout
