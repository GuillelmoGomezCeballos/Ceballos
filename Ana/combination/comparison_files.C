void comparison_files(){
  TFile *_file0 = TFile::Open("/afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/alphaCards/125/hww-12.10fb.mH125.of_1j_shape.root");
  _file0->cd();
  TH1D * h0 = histo_WW->Clone();
  h0->SetDirectory(0);

  TFile *_file1 = TFile::Open("/afs/cern.ch/user/c/ceballos/releases/CMSSW_4_2_8/src/Smurf/cards/HWW2l/alphaCards/130/hwwof_1j.input_8TeV.root"); 
  _file1->cd();
  TH1D * h1 = histo_qqWW->Clone();
  h1->SetDirectory(0);

  printf("sum: %10.6f - %10.6f = %10.6f\n",h0->GetSumOfWeights(),h1->GetSumOfWeights(),h0->GetSumOfWeights()-h1->GetSumOfWeights());
  double scale = h0->GetSumOfWeights()/h1->GetSumOfWeights();
  for(Int_t i=0;i<=h0->GetNbinsX();i++) {if(TMath::Abs(h0->GetBinContent(i)-h1->GetBinContent(i)) > 0.000001){
    printf("%2d   %10.6f - %10.6f = %10.6f",i,h0->GetBinContent(i),h1->GetBinContent(i)*scale,(h0->GetBinContent(i)-h1->GetBinContent(i)*scale)/h0->GetBinContent(i));
    if(TMath::Abs(h0->GetBinContent(i)-h1->GetBinContent(i)*scale) > 0.000001) printf(" BAD\n"); else printf("\n");}
  }
}
