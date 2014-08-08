#!/bin/bash
#Needs directory to process as a command line argument.
if [ $# -lt 1 ]
then
	echo "Directory to analyse is needed."
	exit 1
fi
INPUTDIR=$1
PROCESS=`echo $INPUTDIR | rev | cut -d '/' -f1 | rev`
echo $PROCESS
QUEUE=short
ATHENADIR=/home/wardrope/athena/17.2.11.4/PhysicsAnalysis/TruthParticleID/McParticleAlgs
OUTPUTDIR=/unix/atlas2/wardrope/2014Paper/HepMCFiles/$PROCESS
ls $OUTPUTDIR
if [ $? -ne 0 ]
then
	echo "$OUTPUTDIR does not exist and will be created."
	mkdir $OUTPUTDIR
fi


ARRAY=(`find $INPUTDIR/*.pool.root`)
if [ ${#ARRAY[*]} -eq 0 ]
then
	echo "Couldn't find files with .pool.root suffix"
fi
for COUNTER in `seq 0 1 ${#ARRAY[*]}`
#for COUNTER in `seq 0 1 0`
do
	echo \#! /bin/sh > $PROCESS$COUNTER.sh
	echo cd /work >> $PROCESS$COUNTER.sh	
	echo mkdir \${PBS_JOBID} >> $PROCESS$COUNTER.sh	
	echo cd \${PBS_JOBID} >> $PROCESS$COUNTER.sh
	#echo export AtlasSetup=/afs/cern.ch/atlas/software/dist/AtlasSetup>> $PROCESS$COUNTER.sh
	#echo alias asetup='source ${AtlasSetup}/scripts/asetup.sh' >> $PROCESS$COUNTER.sh
	echo source /afs/cern.ch/atlas/software/dist/AtlasSetup/scripts/asetup.sh 17.2.11.4,here >> $PROCESS$COUNTER.sh
	echo cp ${ARRAY[$COUNTER]} . >> $PROCESS$COUNTER.sh
	echo THIS=${ARRAY[$COUNTER]} >> $PROCESS$COUNTER.sh
	echo sed s/ssssss/\${THIS##*/}/ $ATHENADIR/share/GenEventAsciiWriter_jobOptions.py \> HepMCWriter_$PROCESS$COUNTER.py >> $PROCESS$COUNTER.sh
	echo athena HepMCWriter_$PROCESS$COUNTER.py >> $PROCESS$COUNTER.sh
	echo mv hepmc.ascii $OUTPUTDIR/${PROCESS}_$COUNTER.hepmc2g >> $PROCESS$COUNTER.sh
	echo cd .. >> $PROCESS$COUNTER.sh
	echo rm -r \${PBS_JOBID} >> $PROCESS$COUNTER.sh
	
	qsub -q $QUEUE $PROCESS$COUNTER.sh
	rm $PROCESS$COUNTER.sh
done
