#! /bin/bash

# set env
source ~/EVAL_SH66 6_1_1;

tag=$1;
mass=$2;
baseDir=$4;

#workspace=$baseDir/$mass/hww_${tag}.text
workspace=$baseDir/$mass/${tag}.text

seed=$RANDOM
mkdir -p /tmp/$USER/${seed}/
cd /tmp/$USER/${seed}/

if [ $3 = "exp" ]; then
    # expected
    echo "combine $workspace -M ProfileLikelihood --expectSignal=1 --significance -t -1 -m ${mass} -n Expected${tag} | tee expected${tag}.txt";
          combine $workspace -M ProfileLikelihood --expectSignal=1 --significance -t -1 -m ${mass} -n Expected${tag} | tee expected${tag}.txt;
    mv higgsCombineExpected${tag}.ProfileLikelihood.mH${mass}.root expected${tag}.txt $baseDir/$mass/;
elif [ $3 = "obs" ]; then
    # observed
    echo "combine $workspace -M ProfileLikelihood --significance -m ${mass} -n Observed${tag} | tee observed${tag}.txt";
          combine $workspace -M ProfileLikelihood --significance -m ${mass} -n Observed${tag} --minimizerAlgo Minuit | tee observed${tag}.txt;
          #combine $workspace -M ProfileLikelihood --significance -m ${mass} -n Observed${tag} --minimizerAlgo Minuit2 -v 9 | tee observed${tag}.txt;
    mv higgsCombineObserved${tag}.ProfileLikelihood.mH${mass}.root observed${tag}.txt $baseDir/$mass/;
else
    echo "don't know what to do for $3";
fi

