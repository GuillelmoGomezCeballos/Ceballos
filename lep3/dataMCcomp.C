#include <iomanip>
#include <iostream>
#include <fstream>
TH1D* Histo;
TH1D* Histomc;

void dataMCcomp(
	       Int_t nsel = 0,
	       Int_t iReBin = 1,
               Char_t myxTitle[] ="",
               Char_t myyTitle[] ="",
	       double xstart = 0.2,
	       Char_t namePlot[] ="", 
               TString fileI ="output_ZHWW.root"
               ) {

bool ApplyNorm = false;
bool MakePRPlot = false;
if(MakePRPlot == true) ApplyNorm = false;

TFile* input = new TFile(fileI);

Int_t option = iReBin;
iReBin = TMath::Abs(iReBin);

char sb[50];
sprintf(sb,"hDVar_0_%d",nsel);
TH1D *HistoS = (TH1D*) gROOT->FindObject(sb);

sprintf(sb,"hDVar_1_%d",nsel);
TH1D *Histo1 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_2_%d",nsel);
TH1D *Histo2 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_3_%d",nsel);
TH1D *Histo3 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_4_%d",nsel);
TH1D *Histo4 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_5_%d",nsel);
TH1D *Histo5 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_6_%d",nsel);
TH1D *Histo6 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_7_%d",nsel);
TH1D *Histo7 = (TH1D*) gROOT->FindObject(sb);

sprintf(sb,"hDVar_8_%d",nsel);
TH1D *HistoH0 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_9_%d",nsel);
TH1D *HistoH1 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_10_%d",nsel);
TH1D *HistoH2 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_11_%d",nsel);
TH1D *HistoH3 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_12_%d",nsel);
TH1D *HistoH4 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_13_%d",nsel);
TH1D *HistoH5 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_14_%d",nsel);
TH1D *HistoH6 = (TH1D*) gROOT->FindObject(sb);
sprintf(sb,"hDVar_15_%d",nsel);
TH1D *HistoH7 = (TH1D*) gROOT->FindObject(sb);

TH1D *HistoB = Histo1->Clone();
HistoB->Add(Histo2);
HistoB->Add(Histo3);
HistoB->Add(Histo4);
HistoB->Add(Histo5);
HistoB->Add(Histo6);
HistoB->Add(Histo7);

TH1D *HistoH = HistoH0->Clone();
HistoH->Add(HistoH1);
HistoH->Add(HistoH2);
HistoH->Add(HistoH3);
HistoH->Add(HistoH4);
HistoH->Add(HistoH5);
HistoH->Add(HistoH6);
HistoH->Add(HistoH7);

printf("Bkgs: %f %f %f %f %f %f %f\n",Histo1->GetSumOfWeights(),Histo2->GetSumOfWeights(),Histo3->GetSumOfWeights(),Histo4->GetSumOfWeights(),Histo5->GetSumOfWeights(),Histo6->GetSumOfWeights(),Histo7->GetSumOfWeights());
HistoS->Rebin(iReBin);
HistoB->Rebin(iReBin);
HistoH->Rebin(iReBin);

Double_t scaleS=HistoS->GetSumOfWeights();
Double_t scaleB=HistoB->GetSumOfWeights();
Double_t scaleH=HistoH->GetSumOfWeights();

printf("Norm S/B: %11.3f + %11.3f / %11.3f = %7.3f -> s/sqrt(b) = %7.3f\n",scaleS,scaleH,scaleB,scaleS/scaleB,scaleS/sqrt(scaleB));

if(ApplyNorm == true){
  HistoS->Scale(1./scaleS);
  HistoB->Scale(1./scaleB);
  HistoH->Scale(1./scaleH);
}

atributes(HistoS,myxTitle,myyTitle,1);
atributes(HistoB,myxTitle,myyTitle,4);
atributes(HistoH,myxTitle,myyTitle,2);

cout << "Black is S, Blue is B, Red is H" << endl;
HistoS->SetDirectory(0);
HistoB->SetDirectory(0);
HistoH->SetDirectory(0);

TCanvas* canvas = new TCanvas("cv", "cv", 500, 500, 500, 500);
canvas->cd();

if(MakePRPlot == false){
  if(option > 0){
    HistoS->Draw("hist");
    HistoB->Draw("same,hist");
    HistoH->Draw("same,hist");
  } else {
    HistoB->Draw("hist");
    HistoS->Draw("same,hist");
    HistoH->Draw("same,hist");
  }
} else {
  HistoB->SetFillStyle(3013);
  HistoH->Add(HistoB);
  HistoS->Add(HistoH);
  HistoS->SetMinimum(0);
  HistoB->SetMinimum(0);
  HistoH->SetMinimum(0);

  HistoS->SetFillStyle(1001);
  HistoS->SetFillColor(0);
  HistoS->Draw("e");
  HistoB->Draw("same,hist");
  HistoS->Draw("same,hist,e");
  HistoB->Draw("same,hist");
}

TLegend* leg = new TLegend(xstart,0.8,xstart+0.2,0.9);							
leg ->SetFillStyle(0);
leg ->SetFillColor(kWhite);
leg ->SetBorderSize(0);
leg->SetTextSize(0.035);
if(MakePRPlot == false) {
  leg->AddEntry(HistoS,"ZH #rightarrow ZWW","f");
  leg->AddEntry(HistoH,"ZH #rightarrow anything else","f");	      
}
else {
  leg->AddEntry(HistoS,"Signal","f");	      
}
leg->AddEntry(HistoB,"All backgrounds","f");
leg->Draw();
labelcms  = new TPaveText(0.60,0.95,0.95,0.98,"NDCBR");
labelcms->SetTextAlign(12);
labelcms->SetTextSize(0.05);
labelcms->SetFillColor(0);
labelcms->AddText("CMS Simulation");
labelcms->SetBorderSize(0);
labelcms->Draw();

canvas->GetFrame()->DrawClone();
canvas->RedrawAxis();
canvas->Update();
if(namePlot != "") canvas->SaveAs(namePlot);

input->Close();

}
// atributes
void atributes(TH1D *histo, Char_t xtitle[]="", Char_t ytitle[]="Fraction", Int_t COLOR = 1){
  histo->ResetAttLine();
  histo->ResetAttFill();
  histo->ResetAttMarker();
  histo->GetYaxis()->SetNdivisions(505);
  histo->GetXaxis()->SetNdivisions(505);
  histo->SetTitle("");
  histo->SetMarkerStyle(21);
  histo->SetMarkerColor(2);
  histo->SetMarkerSize(1.2);
  histo->SetLineWidth(4);
  histo->SetLineColor(COLOR);
  histo->SetMarkerStyle(kFullDotLarge);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->SetTitleOffset(1.);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(kTRUE);
}
