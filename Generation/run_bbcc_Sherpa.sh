#!/bin/bash
QUEUE=medium
SHERPADIR=/unix/atlas3/wardrope/Generation/SHERPA-MC-2.1.1
PROCESSDIR=2307_bbcc_1
OUTPUTDIR=/unix/atlas3/wardrope/2014Paper/HepMCFiles/pp_bbcc/
NUMJOBS=500 #number of jobs
COUNTER=201
while [ $COUNTER -le $NUMJOBS ]
do
	echo \#! /bin/sh > bbcc$COUNTER.sh
	echo cd /work >> bbcc$COUNTER.sh	
	echo mkdir \${PBS_JOBID} >> bbcc$COUNTER.sh	
	echo cd \${PBS_JOBID} >> bbcc$COUNTER.sh
	echo cp $SHERPADIR/$PROCESSDIR/Input.tar.gz . >> bbcc$COUNTER.sh
	echo tar xvzf Input.tar.gz >> bbcc$COUNTER.sh
	echo $SHERPADIR/bin/Sherpa -R \${PBS_JOBID} >> bbcc$COUNTER.sh
	echo cp bbcc.hepmc2g $OUTPUTDIR/bbcc_Sherpa_$COUNTER.hepmc2g >> bbcc$COUNTER.sh
	echo cd .. >> bbcc$COUNTER.sh
	echo rm -r \${PBS_JOBID} >> bbcc$COUNTER.sh
	qsub -q $QUEUE bbcc$COUNTER.sh
	rm bbcc$COUNTER.sh
	let COUNTER=COUNTER+1
done
