#/bin/bash
#Needs directory to process as a command line argument.
if [ $# -lt 1 ]
then
	echo "Directory to analyse is needed."
	exit 1
fi
INPUTDIR=$1
PROCESS=`echo $INPUTDIR | rev | cut -d '/' -f1 | rev`
echo $PROCESS
QUEUE=medium
RIVETDIR=/unix/atlas3/wardrope/2014Paper/Rivet2HHAnalysis
OUTPUTDIR=/unix/atlas3/wardrope/2014Paper/AnaResults/$PROCESS
ls $OUTPUTDIR
if [ $? -ne 0 ]
then
	mkdir OUTPUTDIR
fi

ARRAY=(`find $INPUTDIR/*.hepmc2g`)
if [ ${#ARRAY[*]} -eq 0 ]
then
	echo "Couldn't find files with .hepmc2g suffix"
	ARRAY=(`find $INPUTDIR/*.hepmc`)
	if [ ${#ARRAY[*]} -eq 0 ]
	then
		echo "Couldn't find files with .hepmc suffix either"
		exit 2
	fi
fi
INCR=1
XSEC=0.0 #in pb
if [ $PROCESS == "ttbar" ]
then
	INCR=1
	XSEC=212.07
elif [ $PROCESS == "HH" ] #0.0348 pb * BR(H->bb)^2
then 
	XSEC=0.011586
elif [ $PROCESS == "pp_bbbb" ]
then 
	INCR=4
	XSEC=146.03
elif [ $PROCESS == "pp_bbcc" ]
then 
	XSEC=317.78
elif [ $PROCESS == "Hjj" ] #Hjj: 6.055296 pb * BR(H->bb)
then 
	XSEC=3.4939
elif [ $PROCESS == "Hbb" ] #Hbb: 0.8479 pb * BR(H->bb)
then 
	XSEC=0.48924
elif [ $PROCESS == "ttH" ] #ttH: 2.345e-1 pb * BR(H->bb)
then 
	INCR=1
	XSEC=0.13531
elif [ $PROCESS == "ZH" ] #ZH: 6.152e-2 pb * BR(H->bb) (0.577)
then 
	INCR=5
	XSEC=0.035497
fi
 
#echo "Increment is $INCR"
#for COUNTER in `seq 0 $INCR $INCR`
for COUNTER in `seq 0 $INCR ${#ARRAY[*]}`
do
	echo \#! /bin/sh > $PROCESS$COUNTER.sh
	echo cd /work >> $PROCESS$COUNTER.sh	
	echo mkdir \${PBS_JOBID} >> $PROCESS$COUNTER.sh	
	echo cd \${PBS_JOBID} >> $PROCESS$COUNTER.sh
	echo export RIVET_ANALYSIS_PATH=$RIVETDIR >> $PROCESS$COUNTER.sh
	echo source /usr/local/gcc43/setup.sh >> $PROCESS$COUNTER.sh
	echo export ROOTSYS=/usr/local/root-v5.34.12-gcc4.3/ >> $PROCESS$COUNTER.sh
	echo export LD_LIBRARY_PATH=\${ROOTSYS}/lib:\${LD_LIBRARY_PATH} >> $PROCESS$COUNTER.sh
	echo export PATH=\${ROOTSYS}/bin:\${PATH} >> $PROCESS$COUNTER.sh
	echo export LD_LIBRARY_PATH=/usr/local/python-2.7.6/lib:\${LD_LIBRARY_PATH} >> $PROCESS$COUNTER.sh
	echo export PATH=/usr/local/python-2.7.6/bin/:\${PATH} >> $PROCESS$COUNTER.sh
	echo source /unix/atlas3/wardrope/Generation/rivet/local/rivetenv.sh >> $PROCESS$COUNTER.sh
	echo source /unix/atlas3/wardrope/Generation/rivet/local/yodaenv.sh >> $PROCESS$COUNTER.sh
	RIVETCMD="rivet -a Summer_2014_Study --runname=$PROCESS -H ${PROCESS}_$COUNTER.yoda --cross-section=$XSEC"
	let UP=$COUNTER+$INCR-1
	for IFILE in `seq $COUNTER 1 $UP`
	do
		RIVETCMD="$RIVETCMD ${ARRAY[$IFILE]}"
	done
	echo $RIVETCMD >> $PROCESS$COUNTER.sh
	echo mv ${PROCESS}_$COUNTER.yoda $OUTPUTDIR/. >> $PROCESS$COUNTER.sh
	echo mv MVAResults.root $OUTPUTDIR/${PROCESS}_$COUNTER.root >> $PROCESS$COUNTER.sh
	echo cd .. >> $PROCESS$COUNTER.sh
	echo rm -r \${PBS_JOBID} >> $PROCESS$COUNTER.sh
	qsub -q $QUEUE $PROCESS$COUNTER.sh
	rm $PROCESS$COUNTER.sh
done

