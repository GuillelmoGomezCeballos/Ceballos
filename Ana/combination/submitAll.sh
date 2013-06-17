#/bin/sh

#for chan in 8TeV_ofshape0 8TeV_ofshape1 8TeV_ofshape 8TeV_bdt 8TeV_cut \
#            8TeV_of_cut0 8TeV_of_cut1 8TeV_of_cut2 8TeV_sf_cut0 8TeV_sf_cut1 8TeV_sf_cut2 \
#            7TeV_ofshape0 7TeV_ofshape1 7TeV_ofshape 7TeV_cut 7TeV_bdt \
#	     7TeV_of_cut0 7TeV_of_cut1 7TeV_of_cut2 7TeV_sf_cut0 7TeV_sf_cut1 7TeV_sf_cut2 \
#	     7p8TeV 7p8TeV_cut
for chan in 8TeV_bdt 8TeV_cut \
            7TeV_cut 7TeV_bdt \
	    7p8TeV 7p8TeV_cut
do
./submit.sh /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/$1 $chan exp;
./submit.sh /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/$1 $chan obs;
./submit.sh /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/$1 $chan lim;
#./submit.sh /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/$1 $chan mlf;
done
