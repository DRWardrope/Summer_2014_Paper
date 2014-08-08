#!/bin/bash
scl enable devtoolset-2 bash
export LD_LIBRARY_PATH=/usr/local/python-2.7.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/boost-1.54.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/unix/atlas3/wardrope/Generation/LHAPDF/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/unix/atlas3/wardrope/Generation/HepMC/lib:$LD_LIBRARY_PATH
export PYTHIA8DATA=/unix/atlas3/wardrope/Generation/pythia8185/xmldoc
QUEUE=medium
PWHGDIR=/unix/atlas3/wardrope/Generation/POWHEG-BOX/hvq/LHC_14TeV
PYTHIADIR=/unix/atlas3/wardrope/Generation/pythia8185/examples
#PROCESSDIR=2207_tt_1
OUTPUTDIR=/unix/atlas1/wardrope/2014Paper/HepMCFiles/ttbar
NUMJOBS=242 #number of jobs
COUNTER=101

while [ $COUNTER -le $NUMJOBS ]
do
	SUFFIX=$COUNTER
	if [ $COUNTER -lt 10 ]
	then
		SUFFIX=000$COUNTER
	elif [ $COUNTER -lt 100 ]
	then
		SUFFIX=00$COUNTER
	elif [ $COUNTER -lt 1000 ]
	then
		SUFFIX=0$COUNTER
fi
#echo $COUNTER "==>" $SUFFIX
#POWHEG BIT
	echo \#! /bin/sh > tt$SUFFIX.sh
    echo cd /work >> tt$SUFFIX.sh	
	echo mkdir \${PBS_JOBID} >> tt$SUFFIX.sh	
	echo cd \${PBS_JOBID} >> tt$SUFFIX.sh
	echo cp $PWHGDIR/pwgseeds.dat $PWHGDIR/powheg.input $PWHGDIR/pwgcounters$SUFFIX.dat $PWHGDIR/pwggrid-$SUFFIX.dat $PWHGDIR/pwgstat-$SUFFIX.dat $PWHGDIR/pwgxgrid.dat . >> tt$SUFFIX.sh
	printf "echo $SUFFIX | $PWHGDIR/../pwhg_main\n" >> tt$SUFFIX.sh
##PYTHIA BIT
	echo ls >> tt$SUFFIX.sh
	echo echo \${LD_LIBRARY_PATH} >> tt$SUFFIX.sh
	echo sed s/ssssss/pwgevents-$SUFFIX.lhe/ $PYTHIADIR/fastjetFilter.cmnd \> filter.cmnd >> tt$SUFFIX.sh
	echo $PYTHIADIR/powhegFilterHepMC.exe filter.cmnd $OUTPUTDIR/tt$SUFFIX.hepmc >> tt$SUFFIX.sh
	echo cd .. >> tt$SUFFIX.sh
	echo rm -r \${PBS_JOBID} >> tt$SUFFIX.sh
	qsub -V -q $QUEUE tt$SUFFIX.sh
	rm tt$SUFFIX.sh
	let COUNTER=COUNTER+1
done
