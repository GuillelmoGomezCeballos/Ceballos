void fake(TString input = "histo_wjets_madgraph_all.root"){
TFile *_file0 = TFile::Open(input);
_file0->cd();

double fake0[4] = {hDWlnFakeSel_4->GetSumOfWeights(),hDWlnFakeSel_5->GetSumOfWeights(),
                   hDWlnFakeSel_6->GetSumOfWeights(),hDWlnFakeSel_7->GetSumOfWeights()};
double fake1[4] = {hDWlnFakeSel_104->GetSumOfWeights(),hDWlnFakeSel_105->GetSumOfWeights(),
                   hDWlnFakeSel_106->GetSumOfWeights(),hDWlnFakeSel_107->GetSumOfWeights()};
double fake2[4] = {hDWlnFakeSel_204->GetSumOfWeights(),hDWlnFakeSel_205->GetSumOfWeights(),
                   hDWlnFakeSel_206->GetSumOfWeights(),hDWlnFakeSel_207->GetSumOfWeights()};
double fake3[4] = {hDWlnFakeSel_304->GetSumOfWeights(),hDWlnFakeSel_305->GetSumOfWeights(),
                   hDWlnFakeSel_306->GetSumOfWeights(),hDWlnFakeSel_307->GetSumOfWeights()};

printf("(f0) W+f+:%9.1f W-f-:%9.1f W+f-:%9.1f W-f+:%9.1f, SS:%9.1f OS:%9.1f SS+OS:%9.1f, W+f+/W+f-:%6.3f W-f-/W-f+:%6.3f SS/OS:%6.3f\n",fake0[0],fake0[1],fake0[2],fake0[3],fake0[0]+fake0[1],fake0[2]+fake0[3],fake0[0]+fake0[1]+fake0[2]+fake0[3],fake0[0]/fake0[2],fake0[1]/fake0[3],(fake0[0]+fake0[1])/(fake0[2]+fake0[3]));
printf("(f1) W+f+:%9.1f W-f-:%9.1f W+f-:%9.1f W-f+:%9.1f, SS:%9.1f OS:%9.1f SS+OS:%9.1f, W+f+/W+f-:%6.3f W-f-/W-f+:%6.3f SS/OS:%6.3f\n",fake1[0],fake1[1],fake1[2],fake1[3],fake1[0]+fake1[1],fake1[2]+fake1[3],fake1[0]+fake1[1]+fake1[2]+fake1[3],fake1[0]/fake1[2],fake1[1]/fake1[3],(fake1[0]+fake1[1])/(fake1[2]+fake1[3]));
printf("(f2) W+f+:%9.1f W-f-:%9.1f W+f-:%9.1f W-f+:%9.1f, SS:%9.1f OS:%9.1f SS+OS:%9.1f, W+f+/W+f-:%6.3f W-f-/W-f+:%6.3f SS/OS:%6.3f\n",fake2[0],fake2[1],fake2[2],fake2[3],fake2[0]+fake2[1],fake2[2]+fake2[3],fake2[0]+fake2[1]+fake2[2]+fake2[3],fake2[0]/fake2[2],fake2[1]/fake2[3],(fake2[0]+fake2[1])/(fake2[2]+fake2[3]));
printf("(f3) W+f+:%9.1f W-f-:%9.1f W+f-:%9.1f W-f+:%9.1f, SS:%9.1f OS:%9.1f SS+OS:%9.1f, W+f+/W+f-:%6.3f W-f-/W-f+:%6.3f SS/OS:%6.3f\n",fake3[0],fake3[1],fake3[2],fake3[3],fake3[0]+fake3[1],fake3[2]+fake3[3],fake3[0]+fake3[1]+fake3[2]+fake2[3],fake3[0]/fake3[2],fake3[1]/fake3[3],(fake3[0]+fake3[1])/(fake3[2]+fake3[3]));
}
