#!/bin/bash
#Needs directory to process as a command line argument.
if [ $# -lt 1 ]
then
	echo "Directory to analyse is needed."
	exit 1
fi
INPUTDIR=$1
QUEUE=medium
RIVETDIR=/unix/atlas3/wardrope/2014Paper/Summer_2014_Paper/RivetAnalysis
OUTPUTDIR=/unix/atlas3/wardrope/2014Paper/AnaResults/BoostedVariationsOnlyBJets
ls $OUTPUTDIR
if [ $? -ne 0 ]
then
	echo "$OUTPUTDIR does not exist and will be created."
	mkdir $OUTPUTDIR
fi


ARRAY=(`find $INPUTDIR/*.hepmc2g`)
if [ ${#ARRAY[*]} -eq 0 ]
then
	echo "Couldn't find files with .hepmc2g suffix"
fi
for COUNTER in `seq 0 1 ${#ARRAY[*]}`
#for COUNTER in `seq 0 1 0`
do
	FILENAME=`echo ${ARRAY[$COUNTER]} | rev | cut -d '/' -f1 | cut -d '.' -f2 | rev`
	echo $FILENAME
	echo \#! /bin/sh > $FILENAME.sh
	echo cd /work >> $FILENAME.sh	
	echo mkdir \${PBS_JOBID} >> $FILENAME.sh	
	echo cd \${PBS_JOBID} >> $FILENAME.sh
	#echo cp ${ARRAY[$COUNTER]} . >> $FILENAME.sh
	echo export RIVET_ANALYSIS_PATH=$RIVETDIR >> $FILENAME.sh
	echo source /usr/local/gcc43/setup.sh >> $FILENAME.sh
	echo export ROOTSYS=/usr/local/root-v5.34.12-gcc4.3/ >> $FILENAME.sh
	echo export LD_LIBRARY_PATH=\${ROOTSYS}/lib:\${LD_LIBRARY_PATH} >> $FILENAME.sh
	echo export PATH=\${ROOTSYS}/bin:\${PATH} >> $FILENAME.sh
	echo export LD_LIBRARY_PATH=/usr/local/python-2.7.6/lib:\${LD_LIBRARY_PATH} >> $FILENAME.sh
	echo export PATH=/usr/local/python-2.7.6/bin/:\${PATH} >> $FILENAME.sh
	echo source /unix/atlas3/wardrope/Generation/rivet/local/rivetenv.sh >> $FILENAME.sh
	echo source /unix/atlas3/wardrope/Generation/rivet/local/yodaenv.sh >> $FILENAME.sh
	#echo rivet -a Summer_2014_Study --runname=$FILENAME -H ${PROCESS}_$COUNTER.yoda ${ARRAY[$COUNTER]} >> $FILENAME.sh
	echo rivet -a Boosted_2012_Analysis --runname=$FILENAME -H ${PROCESS}_$COUNTER.yoda ${ARRAY[$COUNTER]} >> $FILENAME.sh
	echo mv ${PROCESS}_$COUNTER.yoda $OUTPUTDIR/$FILENAME.yoda >> $FILENAME.sh
	echo mv MVAResults.root $OUTPUTDIR/$FILENAME.root >> $FILENAME.sh
	echo cd .. >> $FILENAME.sh
	echo rm -r \${PBS_JOBID} >> $FILENAME.sh
	qsub -q $QUEUE $FILENAME.sh
	rm $FILENAME.sh
done
