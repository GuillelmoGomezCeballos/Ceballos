#! /bin/bash

LogDir=lsf;
mkdir -p $LogDir;

baseDir=$1
ls $baseDir >& /dev/null;
if [ $? != 0 ]; then
  echo "Wrong folder";
  echo "e.g. /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/ana_test";
  exit;
fi

tag=$2

for mass in 110 115 120 125 130 135 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  if [ $3 = exp ] || [ $3 = obs ] 
      then
      echo "bsub -q cmscaf1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/${3}_${mass}${tag}_${tag}.err computeSignificance.sh ${tag} ${mass} ${3} ${baseDir}"
            bsub -q       1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/${3}_${mass}${tag}_${tag}.err computeSignificance.sh ${tag} ${mass} ${3} ${baseDir}
  elif [ $3 = mlf ] 
      then
      echo "bsub -q cmscaf1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/$3_${mass}_${tag}.err computeMuHatHWW.sh ${tag} ${mass} ${baseDir}"
            bsub -q       1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/$3_${mass}_${tag}.err computeMuHatHWW.sh ${tag} ${mass} ${baseDir}
  elif [ $3 = lim ] 
      then
      echo "bsub -q cmscaf1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/$3_${mass}_${tag}.err computeLimit.sh ${tag} ${mass} ${baseDir}"
            bsub -q       1nd -o ${LogDir}/${3}_${mass}_${tag}.out -e ${LogDir}/$3_${mass}_${tag}.err computeLimit.sh ${tag} ${mass} ${baseDir}
  else exit
  fi
done
