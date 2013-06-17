#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"


TString format = "pdf";


enum {right, top};


//------------------------------------------------------------------------------
// HiggsPlots
//------------------------------------------------------------------------------
void HiggsPlots()
{
  gInterpreter->ExecuteMacro("HiggsPaperStyle.C");

  gSystem->mkdir(format, kTRUE);

  DrawHistogram("wwpresel_0j_mh125_massll",	"m_{#font[12]{ll}}",	      10, 0, "GeV",   right, 125, "0-jet");
  DrawHistogram("wwpresel_1j_mh125_massll",	"m_{#font[12]{ll}}",	      10, 0, "GeV",   right, 125, "1-jet");
  DrawHistogram("wwpresel_0j_mh125_deltaphill", "#Delta#phi_{#font[12]{ll}}", 10, 0, "#circ", top  , 125, "0-jet");
  DrawHistogram("wwpresel_1j_mh125_deltaphill", "#Delta#phi_{#font[12]{ll}}", 10, 0, "#circ", top  , 125, "1-jet");
  DrawHistogram("hwwsel_0j_mh125_massem",	"m_{e#mu}",		      10, 0, "GeV",   right, 125, "0-jet");
  DrawHistogram("hwwsel_1j_mh125_massem",	"m_{e#mu}",		      10, 0, "GeV",   right, 125, "1-jet");

  DrawHistogram("hww_highmass_mll500_0j",     "m_{#font[12]{ll}}",          10, 0, "GeV",   right, 500, "0-jet", 1);
  DrawHistogram("hww_highmass_mll500_1j",     "m_{#font[12]{ll}}",          10, 0, "GeV",   right, 500, "0-jet", 1);
  DrawHistogram("hww_highmass_bdt500_0j",     "BDT Output",	            10, 1, "NULL",  right, 500, "0-jet", 1);
  DrawHistogram("hww_highmass_bdt500_1j",     "BDT Output",	            10, 1, "NULL",  right, 500, "0-jet", 1);

}


//------------------------------------------------------------------------------
// DrawHistogram
//------------------------------------------------------------------------------
void DrawHistogram(TString  input,
		   TString  xtitle         = "xtitle",
		   Int_t    ngroup         = -1,
		   Int_t    precision      = 1,
		   TString  units          = "NULL",
		   Int_t    legendPosition = right,
		   Int_t    mH             = 125,
		   TString  commentData    = "",
		   Bool_t   isLogY         = false,
		   Double_t xmin           = -999,
		   Double_t xmax           =  999,
		   Bool_t   moveOverflow   = true)
{		   
  // Read input file
  //----------------------------------------------------------------------------
  TFile* file = new TFile(Form("inputs/%s.root", input.Data()));

  TH1F* hWW     = (TH1F*)file->Get("histo0");
  TH1F* hZj     = (TH1F*)file->Get("histo1");
  TH1F* htt     = (TH1F*)file->Get("histo2");
  TH1F* hWZ     = (TH1F*)file->Get("histo3");
  TH1F* hWj     = (TH1F*)file->Get("histo4");
  TH1F* hData   = (TH1F*)file->Get("histo5");
  TH1F* hSignal = (TH1F*)file->Get("histos");

  
  // Rebin
  //----------------------------------------------------------------------------
  if (ngroup > 0) {
    hWW    ->Rebin(ngroup);
    hZj    ->Rebin(ngroup);
    htt    ->Rebin(ngroup);
    hWZ    ->Rebin(ngroup);
    hWj    ->Rebin(ngroup);
    hData  ->Rebin(ngroup);
    hSignal->Rebin(ngroup);
  }


  // Show underflow and overflow bins
  //----------------------------------------------------------------------------
  if (moveOverflow) {
    MoveOverflowBins(hWW,     xmin, xmax);
    MoveOverflowBins(hZj,     xmin, xmax);
    MoveOverflowBins(htt,     xmin, xmax);
    MoveOverflowBins(hWZ,     xmin, xmax);
    MoveOverflowBins(hWj,     xmin, xmax);
    MoveOverflowBins(hData,   xmin, xmax);
    MoveOverflowBins(hSignal, xmin, xmax);
  }


  // Stack
  //----------------------------------------------------------------------------
  THStack* hs = new THStack();
  
  hs->Add(hWW);
  hs->Add(hZj);
  hs->Add(htt);
  hs->Add(hWZ);
  hs->Add(hWj);
  hs->Add(hSignal);


  // Cosmetics
  //----------------------------------------------------------------------------
  TH1Cosmetics(hWW,     color_WW,    1001, color_WW,    1);
  TH1Cosmetics(hZj,     color_Zj,    1001, color_Zj,    1);
  TH1Cosmetics(htt,     color_tt,    1001, color_tt,    1);
  TH1Cosmetics(hWZ,     color_WZ,    1001, color_WZ,    1);
  TH1Cosmetics(hWj,     color_Wj,    1001, color_Wj,    1);
  TH1Cosmetics(hSignal, color_higgs, 3354, color_higgs, 1);

  hData->SetLineWidth (1);
  hData->SetMarkerSize(1.1);
  hData->SetLineColor  (kBlack);
  hData->SetMarkerStyle(kFullCircle);

  // Draw
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas(input, input);
  if(isLogY == true) canvas->SetLogy();

  hData->Draw("ep");
  hs   ->Draw("hist,same");
  hData->Draw("ep,same");


  // Axis labels
  //----------------------------------------------------------------------------
  TAxis* xaxis = hData->GetXaxis();
  TAxis* yaxis = hData->GetYaxis();
  
  xaxis->SetNdivisions(505);
  yaxis->SetNdivisions(505);

  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitleOffset(1.5);

  TString ytitle = Form("Entries / %s.%df", "%", precision);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(Form(ytitle.Data(), hData->GetBinWidth(0)));

  if (!units.Contains("NULL")) {

    xaxis->SetTitle(Form("%s (%s)", xaxis->GetTitle(), units.Data()));

    if (units.Contains("circ"))
      yaxis->SetTitle(Form("%s%s", yaxis->GetTitle(), units.Data()));
    else
      yaxis->SetTitle(Form("%s %s", yaxis->GetTitle(), units.Data()));
  }


  // Adjust scale
  //----------------------------------------------------------------------------
  Float_t theMax   = GetMaximumIncludingErrors(hData);
  Float_t theMaxMC = GetMaximumIncludingErrors((TH1F*)hs->GetHistogram());

  if (theMaxMC > theMax) theMax = theMaxMC;

  if (canvas->GetLogy()) {
    hData->SetMaximum(500 * theMax);
    hData->SetMinimum(0.05);
  } else {
    hData->SetMaximum(1.5 * theMax);
  }

  canvas->Modified();
  canvas->Update();


  // Legend
  //----------------------------------------------------------------------------
  Double_t yoffset = 0.048;

  Double_t x0 = 0.610; 
  Double_t y0 = 0.835;
  Double_t x1 = x0; 
  Double_t y1 = y0;
  Double_t x2 = x0; 
  Double_t y2 = y0;

  if (legendPosition == top) {

    x0 -= 0.375;
    x1 -= 0.055;
    x2 += 0.100;

    y1 += 3.*(yoffset+0.001);
    y2 += 5.*(yoffset+0.001);
  }

  DrawTLegend(x0, y0,                      hData,   Form("%s data",commentData.Data()), "ep", 0.035, 0.2, yoffset);
  DrawTLegend(x0, y0 - 1.*(yoffset+0.001), hSignal, Form("m_{H} = %d GeV",mH), "f",  0.035, 0.2, yoffset);
  DrawTLegend(x0, y0 - 2.*(yoffset+0.001), hWW,     " WW",              "f",  0.035, 0.2, yoffset);
  DrawTLegend(x1, y1 - 3.*(yoffset+0.001), hWZ,     " VV",              "f",  0.035, 0.2, yoffset);
  DrawTLegend(x1, y1 - 4.*(yoffset+0.001), htt,     " top",             "f",  0.035, 0.2, yoffset);
  DrawTLegend(x2, y2 - 5.*(yoffset+0.001), hZj,     " Z+jets",          "f",  0.035, 0.2, yoffset);
  DrawTLegend(x2, y2 - 6.*(yoffset+0.001), hWj,     " W+jets",          "f",  0.035, 0.2, yoffset);


  // CMS titles
  //----------------------------------------------------------------------------
  hData->SetTitle("");

  DrawTLatex(0.185, 0.975, 0.05, 13, "CMS");
  DrawTLatex(0.940, 0.983, 0.05, 33, "#sqrt{s} = 8 TeV, L = 5.1 fb^{-1}");


  // Save
  //----------------------------------------------------------------------------
  canvas->GetFrame()->DrawClone();
  canvas->RedrawAxis();
  canvas->Update();
  
  canvas->SaveAs(Form("%s/%s.%s", format.Data(), input.Data(), format.Data()));
}


//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F* h)
{
  Float_t maxWithErrors = 0;

  for (Int_t i=1; i<=h->GetNbinsX(); i++) {

    Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

    if (binHeight > maxWithErrors) maxWithErrors = binHeight;
  }

  return maxWithErrors;
}


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t    x,
		Double_t    y,
		Double_t    tsize,
		Short_t     align,
		const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC();
  tl->SetTextAlign(align);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}


//------------------------------------------------------------------------------
// DrawTLegend
//------------------------------------------------------------------------------
TLegend* DrawTLegend(Float_t x1,
		     Float_t y1,
		     TH1*    hist,
		     TString label,
		     TString option,
		     Float_t tsize   = 0.04,
		     Float_t xoffset = 0.25,
		     Float_t yoffset = 0.07)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + xoffset,
				y1 + yoffset);
  
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();

  return legend;
}


//------------------------------------------------------------------------------
// TH1Cosmetics
//------------------------------------------------------------------------------
void TH1Cosmetics(TH1*    hist,
		  Color_t fcolor,
		  Style_t fstyle,
		  Color_t lcolor,
		  Width_t lwidth)
{
  hist->SetFillColor(fcolor);
  hist->SetFillStyle(fstyle);
  hist->SetLineColor(lcolor);
  hist->SetLineWidth(lwidth);   
}


//------------------------------------------------------------------------------
// MoveOverflowBins
//------------------------------------------------------------------------------
void MoveOverflowBins(TH1* h, Double_t xmin, Double_t xmax) const
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin > -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax <  999) ? axis->FindBin(xmax) : nbins;

  Double_t firstVal = 0;
  Double_t firstErr = 0;

  Double_t lastVal = 0;
  Double_t lastErr = 0;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i <= firstBin) {
      firstVal += h->GetBinContent(i);
      firstErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i >= lastBin) {
      lastVal += h->GetBinContent(i);
      lastErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }

  firstErr = sqrt(firstErr);
  lastErr  = sqrt(lastErr);

  h->SetBinContent(firstBin, firstVal);
  h->SetBinError  (firstBin, firstErr);

  h->SetBinContent(lastBin, lastVal);
  h->SetBinError  (lastBin, lastErr);
}
