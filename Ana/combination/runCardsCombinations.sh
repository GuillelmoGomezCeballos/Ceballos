#! /bin/bash

source ~/EVAL_SH66 6_1_1;

baseDir=$1;
ls $baseDir >& /dev/null;
if [ $? != 0 ]; then
  echo "Wrong folder";
  echo "e.g. /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/ana_test";
  exit;
fi

for mass in 110 115 120 125 130 135 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
  do
  cd $baseDir/${mass}
  echo "combine cards for MH=${mass} "
  # spin
  # combineCards.py -S hwwof_0j_8TeV.txt hwwof_1j_8TeV.txt hwwof_0j_7TeV.txt hwwof_1j_7TeV.txt > hww_spin.text
  # 8 TeV
  combineCards.py -S hwwof_0j_shape_8TeV.txt > hww_8TeV_ofshape0.text
  combineCards.py -S hwwof_1j_shape_8TeV.txt > hww_8TeV_ofshape1.text
  combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt  >  hww_8TeV_ofshape.text
  combineCards.py -S hwwof_0j_cut_8TeV.txt  >  hww_8TeV_of_cut0.text
  combineCards.py -S hwwof_1j_cut_8TeV.txt  >  hww_8TeV_of_cut1.text
  combineCards.py -S hwwof_2j_cut_8TeV.txt  >  hww_8TeV_of_cut2.text
  combineCards.py -S hwwsf_0j_cut_8TeV.txt  >  hww_8TeV_sf_cut0.text
  combineCards.py -S hwwsf_1j_cut_8TeV.txt  >  hww_8TeV_sf_cut1.text
  combineCards.py -S hwwsf_2j_cut_8TeV.txt  >  hww_8TeV_sf_cut2.text
  #combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt >  hww_8TeV_bdt.text
  combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt >  hww_8TeV_bdt.text
  combineCards.py -S hwwof_0j_cut_8TeV.txt   hwwof_1j_cut_8TeV.txt   hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt >  hww_8TeV_cut.text
  # 7 TeV
  combineCards.py -S hwwof_0j_shape_7TeV.txt > hww_7TeV_ofshape0.text
  combineCards.py -S hwwof_1j_shape_7TeV.txt > hww_7TeV_ofshape1.text
  combineCards.py -S hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt  >  hww_7TeV_ofshape.text
  combineCards.py -S hwwof_0j_cut_7TeV.txt  >  hww_7TeV_of_cut0.text
  combineCards.py -S hwwof_1j_cut_7TeV.txt  >  hww_7TeV_of_cut1.text
  combineCards.py -S hwwof_2j_cut_7TeV.txt  >  hww_7TeV_of_cut2.text
  combineCards.py -S hwwsf_0j_cut_7TeV.txt  >  hww_7TeV_sf_cut0.text
  combineCards.py -S hwwsf_1j_cut_7TeV.txt  >  hww_7TeV_sf_cut1.text
  combineCards.py -S hwwsf_2j_cut_7TeV.txt  >  hww_7TeV_sf_cut2.text
  #combineCards.py -S hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt hwwsf_2j_cut_7TeV.txt hwwof_2j_cut_7TeV.txt >  hww_7TeV_bdt.text
  combineCards.py -S hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_7TeV_bdt.text
  combineCards.py -S hwwof_0j_cut_7TeV.txt   hwwof_1j_cut_7TeV.txt   hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_7TeV_cut.text
  # 7+8 TeV
  #combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt hwwof_2j_cut_8TeV.txt hwwsf_2j_cut_8TeV.txt \
  #                   hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt hwwsf_2j_cut_7TeV.txt hwwof_2j_cut_7TeV.txt >  hww_7p8TeV.text
  combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt \
                     hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_7p8TeV.text
  combineCards.py -S hwwof_0j_cut_8TeV.txt   hwwof_1j_cut_8TeV.txt   hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt \
                     hwwof_0j_cut_7TeV.txt   hwwof_1j_cut_7TeV.txt   hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_7p8TeV_cut.text
  # VBF analysis
  combineCards.py -S hwwsf_2j_cut_?TeV.txt hwwof_2j_cut_?TeV.txt   > hww_7p8TeV_vbf_cut.text
  combineCards.py -S hwwsf_2j_cut_?TeV.txt hwwof_2j_shape_?TeV.txt > hww_7p8TeV_vbf_bdt.text
  # WH->3l
  #combineCards.py -S vh3l1_cut_7TeV.txt   vh3l2_cut_7TeV.txt   vh3l1_cut_8TeV.txt   vh3l2_cut_8TeV.txt   > hww_wh_3l_cut.text
  #combineCards.py -S vh3l1_shape_7TeV.txt vh3l2_shape_7TeV.txt vh3l1_shape_8TeV.txt vh3l2_shape_8TeV.txt > hww_wh_3l_bdt.text
done
