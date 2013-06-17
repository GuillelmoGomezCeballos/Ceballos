#! /bin/bash

#combine -M AsymptoticNew ${workspace} -m ${mass} --scanMin=0.01 --scanMax=5 --nPoints=100 -n hWW${tag} | tee limit${tag}.txt;

# set env
source ~/EVAL_SH66 6_1_1;

tag=$1;
mass=$2;
baseDir=$3;

#workspace=$baseDir/$mass/hww_${tag}.text;
workspace=$baseDir/$mass/${tag}.text;

seed=$RANDOM
mkdir -p /tmp/$USER/${seed}/;
cd /tmp/$USER/${seed}/;

echo "combine -M Asymptotic ${workspace} --newExpected=1 -m ${mass} -n hWW${tag} | tee limit${tag}.txt";
      combine -M Asymptotic ${workspace} --newExpected=1 -m ${mass} -n hWW${tag} | tee limit${tag}.txt;
mv higgsCombinehWW${tag}.Asymptotic.mH${mass}.root limit${tag}.txt $baseDir/$mass/;
# --expectSignal=1 -t -1
