#!/bin/bash
QUEUE=medium
SHERPADIR=/unix/atlas3/wardrope/Generation/SHERPA-MC-2.1.1
PROCESSDIR=2207_bbbb_1
OUTPUTDIR=/unix/atlas3/wardrope/2014Paper/HepMCFiles/pp_bbbb/
NUMJOBS=700 #number of jobs
COUNTER=501
while [ $COUNTER -le $NUMJOBS ]
do
	echo \#! /bin/sh > bbbb$COUNTER.sh
	echo cd /work >> bbbb$COUNTER.sh	
	echo mkdir \${PBS_JOBID} >> bbbb$COUNTER.sh	
	echo cd \${PBS_JOBID} >> bbbb$COUNTER.sh
	echo cp $SHERPADIR/$PROCESSDIR/Input.tar.gz . >> bbbb$COUNTER.sh
	echo tar xvzf Input.tar.gz >> bbbb$COUNTER.sh
	echo $SHERPADIR/bin/Sherpa -R \${PBS_JOBID} >> bbbb$COUNTER.sh
	echo cp bbbb.hepmc2g $OUTPUTDIR/bbbb_Sherpa_$COUNTER.hepmc2g >> bbbb$COUNTER.sh
	echo cd .. >> bbbb$COUNTER.sh
	echo rm -r \${PBS_JOBID} >> bbbb$COUNTER.sh
	qsub -q $QUEUE bbbb$COUNTER.sh
	rm bbbb$COUNTER.sh
	let COUNTER=COUNTER+1
done
