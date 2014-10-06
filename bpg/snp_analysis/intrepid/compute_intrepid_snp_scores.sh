#!/bin/bash
#PBS -q research
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -t 0-99
#PBS -d /home/curth/snp_analysis/runs/run_intrepid/working_dir
export NUMTHREADS=100
export NUMPROTSPROC=5000
export PYTHONPATH=/usr/lib/python2.4/site-packages:/clusterfs/ohana/software/webserver/test/bin/modeller9v8/lib/x86_64-intel8:/clusterfs/ohana/software/webserver/test/bin/modeller9v8/modlib:/clusterfs/ohana/software/webserver/test/bin/modeller9v8/modlib:/home/curth/ohana_repository:/home/curth/unfuddle/ohana_repository:/home/curth/ohana_repository/bpg/snp_analysis/rostlab:/home/curth/ohana_repository/bpg/snp_analysis/intrepid:/home/curth/ohana_repository/bpg/snp_analysis/common

/home/curth/ohana_repository/bpg/snp_analysis/intrepid/compute_intrepid_snp_scores.py $PBS_ARRAYID $NUMTHREADS $NUMPROTSPROC 
