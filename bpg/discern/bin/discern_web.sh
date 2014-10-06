#!/bin/bash

export LD_LIBRARY_PATH=/clusterfs/ohana/software/BALL-1.2/contrib/lib:/clusterfs/ohana/software/BALL-1.2/lib/Linux-x86_64-g++_4.1.2
export PYTHONPATH=/clusterfs/ohana/software/lib/python2.5/site-packages/BPG_common:/clusterfs/ohana/software/BALL-1.2//contrib/lib:/clusterfs/ohana/software/BALL-1.2//lib/Linux-x86_64-g++_4.1.2
export BALL=/clusterfs/ohana/software/BALL-1.2/

BINDIR=/clusterfs/ohana/software/discern/bin

${BINDIR}/pipeline_web.py $1 ${BINDIR}/parameters.train-all/weights.txt ${BINDIR}/parameters.train-all/mean.txt ${BINDIR}/parameters.train-all/stddev.txt

