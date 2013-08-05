#/bin/sh

export SEL=$1;
export NJETS=$2;

if [ $SEL == 3 ] && [ $NJETS == 0 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh105inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,0,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh115inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,0,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh125inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,0,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh135inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,0,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh145inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,0,120,160,0.25,180)';

elif [ $SEL == 3 ] && [ $NJETS == 1 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh105inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,1,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh115inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,1,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh125inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,1,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh135inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,1,120,160,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_53x/backgroundA_skim10.root","ntuples_53x/zh145inv.root","ntuples_53x/data_skim10.root","ntuples_53x/hww_syst_skim10_gj.root",3,2,1,120,160,0.25,180)';

elif [ $SEL == 4 ] && [ $NJETS == 0 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh105inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,0,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh115inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,0,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh125inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,0,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh135inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,0,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh145inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,0,120,160,0.25,180)'

elif [ $SEL == 4 ] && [ $NJETS == 1 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh105inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,1,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh115inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,1,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh125inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,1,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh135inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,1,120,160,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_42x_v9/backgroundZH_skim10.root","ntuples_42x_v9/zh145inv.root","ntuples_42x_v9/data_skim10.root","ntuples_42x_v9/zh_syst_skim10_gj.root",4,2,1,120,160,0.25,180)'

fi
