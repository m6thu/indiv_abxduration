#!/bin/sh
#
# EXAMPLE MPICH SCRIPT FOR SGE
#
# Name of the job
#$ -N test_opencv
#
# pe request for MPICH. Set your number of processors here.

# Make sure you use the "mpich" parallel environment.
#$ -pe mpich 1
#
# Set to working directory to current directory
#$ -cwd
#
# Combine stderr and stdout to one file
#$ -j y
#
# Run job through bash shell
#$ -S /bin/bash
#
# Adjust MPICH procgroup to ensure smooth shutdown
export MPICH_PROCESS_GROUP=no
#
# Print useful information
echo "Program start at: `/bin/date`"
echo "Got $NSLOTS slots."
echo "Machines:"
cat $TMPDIR/machines

# Run program. Should use full path instead of relative path
#/usr/bin/R --vanilla < /home/Mathupanee/grid_model_cluster/tests/script_xcorrelate_test.R
echo /home/Mathupanee/grid_model_cluster/run_cluster_simple.sh
echo "Program finish at: ‘/bin/date‘"
