#!/bin/tcsh

# Read the arguments
echo "Starting catalogFiles.csh with arguments:  $*"
set isAll=$1
set catalogDir=$2
set book=$3
set dataset=$4
set channel=$6
set events=$7
set noStage=$8
set version=1

# Prepare environment
echo " "
echo "isAll     : " $isAll
echo "catalogDir: " $catalogDir
echo "book      : " $book
echo "dataset   : " $dataset
echo "channel   : " $channel
echo "events    : " $events
echo "noStage   : " $noStage
echo "version   : " $version
echo " "

set STAGE_SVCCLASS=cmsprod

set workDir=$HOME/condor
cd $workDir
set script=run.csh

#set dataDir=`tail -1  $catalogDir/$book/$dataset/Filesets | cut -d' ' -f2`
set dataDir=`tail -1  $catalogDir/$book/$dataset/Filesets | awk '{print$2}'`

# Stage the entire dataset to be processed
if ( "$noStage" == "" ) then
  set var=`echo $dataDir | grep castor/cern.ch`
  if ( "$var" != "" ) then
    echo " "
    echo "Staging files into castor in $dataDir/$book/$dataset"
    echo " "
    /home/mitprod/bin/stageSample.py --dataDir=$dataDir
  else
    echo " "
    echo "No Staging, files are local"
    echo " "
  endif
else
  echo " "
  echo "No pre-staging of the entire sample."
  echo " "
endif

# Looping through each single file and submitting the condor jobs
echo " "
echo "Submitting jobs to condor"
echo " "
echo " cat $catalogDir/$book/$dataset/Filesets | grep -v ^# | cut -d' ' -f1 "

set skimDir=""

set logFile=/tmp/file_${skimDir}${channel}_${events}_${version}.log
rm -f submit.cmd

set runMacro="runHWWSelection.C+"
if ( "${skimDir}" != "" ) then
  set runMacro="doubleLepton.C+"
  #set runMacro="TightPlusDenominatorSkim.C+"
  mkdir -p ${skimDir}/$dataset;
endif

#foreach file (`cat $catalogDir/$book/$dataset/Filesets | cut -d' ' -f1 `)
foreach file  ($5)
  #echo $file
  set logFile=`echo $book/$dataset/${skimDir}/$file | tr '/' '+'`
  set logFile=/tmp/$logFile
  rm -f $logFile
  echo " $script $catalogDir $book $dataset $file"
  mkdir -p $workDir/test/${skimDir}/${catalogDir}/${book}/${dataset}
  rm    -f $workDir/test/${skimDir}/${catalogDir}/${book}/${dataset}/${file}_${channel}_${events}_${version}.out

  if ( "$isAll" == "1" ) then
    set file = "all"
  endif

  bsub -q 1nd -o $workDir/test/${skimDir}/${catalogDir}/${book}/${dataset}/${file}_${channel}_${events}_${version}.out -J ${dataset}_${file}_${channel} $HOME/condor/run.csh $runMacro $catalogDir $book $dataset $file $channel $events

  if ( "$isAll" == "1" ) then
    exit 0  
  endif
   
end

exit 0
