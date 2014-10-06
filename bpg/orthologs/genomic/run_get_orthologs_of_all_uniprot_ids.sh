#!/bin/bash
#PBS -l nodes=25:ppn=3
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -o run_get_orthologs_of_all_uniprot_ids.out
##PBS -M ruchira@berkeley.edu
##PBS -m abe
#PBS -q library
mpirun -np 75 -hostfile $PBS_NODEFILE -x PYTHONPATH="/home/ruchira/ohana_repository:/home/ruchira/lib64/python" -x LD_LIBRARY_PATH="/usr/lib64/openmpi/1.4-gcc/lib" /home/ruchira/ohana_repository/bpg/orthologs/genomic/get_orthologs_of_all_uniprot_ids.py
