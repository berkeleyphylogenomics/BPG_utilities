#!/bin/bash

#PBS -q library
#PBS -t 0-56
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o run_array_of_write_pairwise_species_file.out

/home/ruchira/ohana_repository/bpg/orthologs/genomic/array_of_write_pairwise_species_file.py
