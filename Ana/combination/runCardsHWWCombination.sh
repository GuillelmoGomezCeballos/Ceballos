#! /bin/bash

source ~/EVAL_SH66 6_1_1;

baseDir=$1;
ls $baseDir >& /dev/null;
if [ $? != 0 ]; then
  echo "Wrong folder";
  echo "e.g. /afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/ana_test";
  exit;
fi

for mass in 110 115 120 125 125.7 130 135 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
  do
  cd $baseDir/${mass}
  echo "combine cards for MH=${mass} "
  combineCards.py -S hwwof_0j_shape_8TeV.txt hwwof_1j_shape_8TeV.txt hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt \
                     hwwof_0j_shape_7TeV.txt hwwof_1j_shape_7TeV.txt hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_01j_shape.text
  combineCards.py -S hwwof_0j_cut_8TeV.txt   hwwof_1j_cut_8TeV.txt   hwwsf_0j_cut_8TeV.txt hwwsf_1j_cut_8TeV.txt \
                     hwwof_0j_cut_7TeV.txt   hwwof_1j_cut_7TeV.txt   hwwsf_0j_cut_7TeV.txt hwwsf_1j_cut_7TeV.txt >  hww_01j_cut.text

  combineCards.py -S hwwsf_2j_cut_?TeV.txt hwwof_2j_shape_?TeV.txt > hww_vbf_shape.text
  combineCards.py -S hwwsf_2j_cut_?TeV.txt hwwof_2j_cut_?TeV.txt   > hww_vbf_cut.text

  combineCards.py -S hwwsf_2j_cut_?TeV.txt   > hww_vbf_sfcut.text
  combineCards.py -S hwwof_2j_cut_?TeV.txt   > hww_vbf_ofcut.text

  if [ ${mass} -lt 210 ]; then
    combineCards.py -S hwwof_?j_shape_?TeV.txt hwwsf_?j_cut_?TeV.txt hww?f_vh2j_cut_?TeV.txt vh3l?_shape_?TeV.txt zh3l2j*shape_???_?TeV.txt > hww_7p8TeV.text

    combineCards.py -S hww?f_vh2j_cut_?TeV.txt > hww_vh2j_cut.text
    # ugly stuff because there are no 7TeV data cards for VH2j in some mass points
    if [ ${mass} -lt 120 ] || [ ${mass} -eq 200 ]; then
      combineCards.py -S hwwsf_vh2j_cut_8TeV.txt hwwof_vh2j_shape_8TeV.txt > hww_vh2j_shape.text
    else
      combineCards.py -S hwwsf_vh2j_cut_?TeV.txt hwwof_vh2j_shape_8TeV.txt hwwof_vh2j_cut_7TeV.txt > hww_vh2j_shape.text
    fi

    combineCards.py -S vh3l?_cut_?TeV.txt   > hww_vh3l_cut.text
    combineCards.py -S vh3l?_shape_?TeV.txt > hww_vh3l_shape.text

    combineCards.py -S zh3l2j*shape_???_?TeV.txt > hww_zh3l2j_shape.text
    combineCards.py -S zh3l2j*cut_???_?TeV.txt   > hww_zh3l2j_cut.text

  elif [ ${mass} -lt 310 ]; then
    combineCards.py -S hwwof_?j_shape_?TeV.txt hwwsf_?j_cut_?TeV.txt hww?f_vh2j_cut_?TeV.txt > hww_7p8TeV.text

    combineCards.py -S hww?f_vh2j_cut_?TeV.txt > hww_vh2j_cut.text
    combineCards.py -S hwwsf_vh2j_cut_8TeV.txt hwwof_vh2j_shape_8TeV.txt > hww_vh2j_shape.text

  else
    combineCards.py -S hwwof_?j_shape_?TeV.txt hwwsf_?j_cut_?TeV.txt > hww_7p8TeV.text
  fi
done
