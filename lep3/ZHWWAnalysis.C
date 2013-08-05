#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TGraphErrors.h>           // Graph class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <set>
#include <vector>
#include <utility>
#include "TH1D.h"
#endif


double DeltaPhi(double phi1, double phi2);
static const Bool_t verbose = kFALSE;

void ZHWWAnalysis(Int_t period = 0)
{

  Double_t lumi = 1.0;
  TString filesPath   = "";
  if     (period == 0){
    lumi = 500;
    filesPath  = "/afs/cern.ch/work/m/mzanetti/public/LEP3/ntuples/WW/";
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  vector<TString> ClassName;
  ClassName.push_back(Form("ZH->qqqqqq"));
  ClassName.push_back(Form("ZH->qqqqln"));
  ClassName.push_back(Form("ZH->qqenmn"));
  ClassName.push_back(Form("ZH->qqlnln"));
  ClassName.push_back(Form("ZH->nnqqqq"));
  ClassName.push_back(Form("ZH->nnenmn"));
  ClassName.push_back(Form("ZH->nnlnln"));
  ClassName.push_back(Form("ZH->llqqqq"));
  ClassName.push_back(Form("ZH->lllnqq"));
  ClassName.push_back(Form("ZH->llenmn"));
  ClassName.push_back(Form("ZH->lllnln"));
  ClassName.push_back(Form("ZH->nnqqln"));
  double nClass[12][10];
  for(int i=0; i<12; i++){
    for(int j=0; j<10; j++){
      nClass[i][j] = 0.0;
    }
  }
  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  // order is critical
  infilenamev.push_back(Form("%s/lep3WWana_signal.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_Zmm.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_Zee.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_Ztt.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_WW.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_ZZ_patrick.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_qqbar.root",filesPath.Data()));
  infilenamev.push_back(Form("%s/lep3WWana_gaga.root",filesPath.Data()));
  //infilenamev.push_back(Form("%s/lep3WWana_eeZ.root",filesPath.Data()));
  //infilenamev.push_back(Form("%s/lep3WWana_enW.root",filesPath.Data()));
  vector<double> scale1fb;  
  scale1fb.push_back(194.0);
  scale1fb.push_back(4200.0);
  scale1fb.push_back(4200.0);
  scale1fb.push_back(4200.0);
  scale1fb.push_back(16000.0);
  scale1fb.push_back(1300.0);
  scale1fb.push_back(50000.0);
  scale1fb.push_back(15000000.0);
  //scale1fb.push_back(3800.0);
  //scale1fb.push_back(1370.0);
  
  
  TH1D *hDSig[10];
  for(UInt_t i=0; i<10; i++){
    hDSig[i] = new TH1D(Form("hDSig_%d" ,i),Form("hDSig_%d" ,i),30,-0.5,29.5);
  }

  TH1D *hDVar[infilenamev.size()+8][100];
  for(UInt_t i=0; i<infilenamev.size()+8; i++){
    hDVar[i][0] = new TH1D(Form("hDVar_%d_0" ,i),Form("hDVar_%d_0" ,i),100,-0.5,99.5);
    hDVar[i][1] = new TH1D(Form("hDVar_%d_1" ,i),Form("hDVar_%d_1" ,i),300,0.,300.);
    hDVar[i][2] = new TH1D(Form("hDVar_%d_2" ,i),Form("hDVar_%d_2" ,i),50,-0.5,49.5);
    hDVar[i][3] = new TH1D(Form("hDVar_%d_3" ,i),Form("hDVar_%d_3" ,i),5,-0.5,4.5);
    hDVar[i][4] = new TH1D(Form("hDVar_%d_4" ,i),Form("hDVar_%d_4" ,i),5,-0.5,4.5);
    hDVar[i][5] = new TH1D(Form("hDVar_%d_5" ,i),Form("hDVar_%d_5" ,i),5,-0.5,4.5);
    hDVar[i][6] = new TH1D(Form("hDVar_%d_6" ,i),Form("hDVar_%d_6" ,i),300,0.,300.);
    hDVar[i][7] = new TH1D(Form("hDVar_%d_7" ,i),Form("hDVar_%d_7" ,i),150,0.,150.);
    hDVar[i][8] = new TH1D(Form("hDVar_%d_8" ,i),Form("hDVar_%d_8" ,i),150,0.,150.);
    hDVar[i][9] = new TH1D(Form("hDVar_%d_9" ,i),Form("hDVar_%d_9" ,i),200,0.,200.);
    hDVar[i][10]= new TH1D(Form("hDVar_%d_10",i),Form("hDVar_%d_10",i),200,0.,200.);
    hDVar[i][11]= new TH1D(Form("hDVar_%d_11",i),Form("hDVar_%d_11",i),200,0.,200.);
    hDVar[i][12]= new TH1D(Form("hDVar_%d_12",i),Form("hDVar_%d_12",i),150,0.,150.);

    hDVar[i][13]= new TH1D(Form("hDVar_%d_13",i),Form("hDVar_%d_13",i),30,0.,3.0);
    hDVar[i][14]= new TH1D(Form("hDVar_%d_14",i),Form("hDVar_%d_14",i),100,-0.5,99.5);
    hDVar[i][15]= new TH1D(Form("hDVar_%d_15",i),Form("hDVar_%d_15",i),300,0.,300.);
    hDVar[i][16]= new TH1D(Form("hDVar_%d_16",i),Form("hDVar_%d_16",i),150,0.,150.);
    hDVar[i][17]= new TH1D(Form("hDVar_%d_17",i),Form("hDVar_%d_17",i),160,0.,160.);
    hDVar[i][18]= new TH1D(Form("hDVar_%d_18",i),Form("hDVar_%d_18",i),200,0.,200.);
    hDVar[i][19]= new TH1D(Form("hDVar_%d_19",i),Form("hDVar_%d_19",i),200,0.,200.);
    hDVar[i][20]= new TH1D(Form("hDVar_%d_20",i),Form("hDVar_%d_20",i),200,0.,200.);
    hDVar[i][21]= new TH1D(Form("hDVar_%d_21",i),Form("hDVar_%d_21",i),200,0.,200.);
    hDVar[i][22]= new TH1D(Form("hDVar_%d_22",i),Form("hDVar_%d_22",i),200,-10.,10.);
    hDVar[i][23]= new TH1D(Form("hDVar_%d_23",i),Form("hDVar_%d_23",i),200,-10.,10.);

    hDVar[i][24]= new TH1D(Form("hDVar_%d_24",i),Form("hDVar_%d_24",i),200,0.,200.);
    hDVar[i][25]= new TH1D(Form("hDVar_%d_25",i),Form("hDVar_%d_25",i),90,0.,180.);
    hDVar[i][26]= new TH1D(Form("hDVar_%d_26",i),Form("hDVar_%d_26",i),90,0.,180.);
    hDVar[i][27]= new TH1D(Form("hDVar_%d_27",i),Form("hDVar_%d_27",i),200,0.,200.);
    hDVar[i][28]= new TH1D(Form("hDVar_%d_28",i),Form("hDVar_%d_28",i),200,0.,200.);
    hDVar[i][29]= new TH1D(Form("hDVar_%d_29",i),Form("hDVar_%d_29",i),200,0.,200.);
    hDVar[i][30]= new TH1D(Form("hDVar_%d_30",i),Form("hDVar_%d_30",i),160,0.,160.);
    hDVar[i][96]= new TH1D(Form("hDVar_%d_96",i),Form("hDVar_%d_96",i),200,-10.,10.);
    hDVar[i][97]= new TH1D(Form("hDVar_%d_97",i),Form("hDVar_%d_97",i),200,-10.,10.);

    hDVar[i][31]= new TH1D(Form("hDVar_%d_31",i),Form("hDVar_%d_31",i),100,-0.5,99.5);
    hDVar[i][32]= new TH1D(Form("hDVar_%d_32",i),Form("hDVar_%d_32",i),100,0.,100.);
    hDVar[i][33]= new TH1D(Form("hDVar_%d_33",i),Form("hDVar_%d_33",i),180,0.,180.);
    hDVar[i][34]= new TH1D(Form("hDVar_%d_34",i),Form("hDVar_%d_34",i),200,0.,200.);
    hDVar[i][35]= new TH1D(Form("hDVar_%d_35",i),Form("hDVar_%d_35",i),200,0.,200.);
    hDVar[i][36]= new TH1D(Form("hDVar_%d_36",i),Form("hDVar_%d_36",i),160,0.,160.);
    hDVar[i][37]= new TH1D(Form("hDVar_%d_37",i),Form("hDVar_%d_37",i),160,0.,160.);

    hDVar[i][38]= new TH1D(Form("hDVar_%d_38",i),Form("hDVar_%d_38",i),200,0.,200.);
    hDVar[i][39]= new TH1D(Form("hDVar_%d_39",i),Form("hDVar_%d_39",i),100,-0.5,99.5);
    hDVar[i][40]= new TH1D(Form("hDVar_%d_40",i),Form("hDVar_%d_40",i),200,0.,200.);
    hDVar[i][41]= new TH1D(Form("hDVar_%d_41",i),Form("hDVar_%d_41",i),100,0.,100.);
    hDVar[i][42]= new TH1D(Form("hDVar_%d_42",i),Form("hDVar_%d_42",i),200,0.,200.);
    hDVar[i][43]= new TH1D(Form("hDVar_%d_43",i),Form("hDVar_%d_43",i),200,0.,200.);
    hDVar[i][44]= new TH1D(Form("hDVar_%d_44",i),Form("hDVar_%d_44",i),200,-10.,10.);
    hDVar[i][45]= new TH1D(Form("hDVar_%d_45",i),Form("hDVar_%d_45",i),200,-10.,10.);
    hDVar[i][46]= new TH1D(Form("hDVar_%d_46",i),Form("hDVar_%d_46",i),160,0.,160.);

    hDVar[i][47]= new TH1D(Form("hDVar_%d_47",i),Form("hDVar_%d_47",i),180,0.,180.);
    hDVar[i][48]= new TH1D(Form("hDVar_%d_48",i),Form("hDVar_%d_48",i),200,0.,200.);
    hDVar[i][49]= new TH1D(Form("hDVar_%d_49",i),Form("hDVar_%d_49",i),200,0.,200.);
    hDVar[i][50]= new TH1D(Form("hDVar_%d_50",i),Form("hDVar_%d_50",i),200,0.,200.);
    hDVar[i][51]= new TH1D(Form("hDVar_%d_51",i),Form("hDVar_%d_51",i),200,0.,200.);
    hDVar[i][52]= new TH1D(Form("hDVar_%d_52",i),Form("hDVar_%d_52",i),180,0.,180.);
    hDVar[i][53]= new TH1D(Form("hDVar_%d_53",i),Form("hDVar_%d_53",i),200,0.,200.);
    hDVar[i][54]= new TH1D(Form("hDVar_%d_54",i),Form("hDVar_%d_54",i),200,0.,200.);
    hDVar[i][55]= new TH1D(Form("hDVar_%d_55",i),Form("hDVar_%d_55",i),200,0.,200.);

    hDVar[i][56]= new TH1D(Form("hDVar_%d_56",i),Form("hDVar_%d_56",i),200,0.,200.);
    hDVar[i][57]= new TH1D(Form("hDVar_%d_57",i),Form("hDVar_%d_57",i),100,-0.5,99.5);
    hDVar[i][58]= new TH1D(Form("hDVar_%d_58",i),Form("hDVar_%d_58",i),200,0.,200.);
    hDVar[i][59]= new TH1D(Form("hDVar_%d_59",i),Form("hDVar_%d_59",i),200,0.,200.);
    hDVar[i][60]= new TH1D(Form("hDVar_%d_60",i),Form("hDVar_%d_60",i),180,0.,180.);
    hDVar[i][61]= new TH1D(Form("hDVar_%d_61",i),Form("hDVar_%d_61",i),200,-10.,10.);
    hDVar[i][62]= new TH1D(Form("hDVar_%d_62",i),Form("hDVar_%d_62",i),200,-10.,10.);
    hDVar[i][63]= new TH1D(Form("hDVar_%d_63",i),Form("hDVar_%d_63",i),160,0.,160.);

    hDVar[i][64]= new TH1D(Form("hDVar_%d_64",i),Form("hDVar_%d_64",i),200,0.,200.);
    hDVar[i][65]= new TH1D(Form("hDVar_%d_65",i),Form("hDVar_%d_65",i),200,0.,200.);
    hDVar[i][66]= new TH1D(Form("hDVar_%d_66",i),Form("hDVar_%d_66",i),180,0.,180.);
    hDVar[i][67]= new TH1D(Form("hDVar_%d_67",i),Form("hDVar_%d_67",i),100,-0.5,99.5);
    hDVar[i][68]= new TH1D(Form("hDVar_%d_68",i),Form("hDVar_%d_68",i),200,0.,200.);
    hDVar[i][69]= new TH1D(Form("hDVar_%d_69",i),Form("hDVar_%d_69",i),200,0.,200.);
    hDVar[i][70]= new TH1D(Form("hDVar_%d_70",i),Form("hDVar_%d_70",i),200,0.,200.);
    hDVar[i][71]= new TH1D(Form("hDVar_%d_71",i),Form("hDVar_%d_71",i),180,0.,180.);
    hDVar[i][72]= new TH1D(Form("hDVar_%d_72",i),Form("hDVar_%d_72",i),160,0.,160.);
    hDVar[i][98]= new TH1D(Form("hDVar_%d_98",i),Form("hDVar_%d_98",i),200,-10.,10.);
    hDVar[i][99]= new TH1D(Form("hDVar_%d_99",i),Form("hDVar_%d_99",i),200,-10.,10.);

    hDVar[i][73]= new TH1D(Form("hDVar_%d_73",i),Form("hDVar_%d_73",i),200,0.,200.);
    hDVar[i][74]= new TH1D(Form("hDVar_%d_74",i),Form("hDVar_%d_74",i),200,0.,200.);
    hDVar[i][75]= new TH1D(Form("hDVar_%d_75",i),Form("hDVar_%d_75",i),200,0.,200.);
    hDVar[i][76]= new TH1D(Form("hDVar_%d_76",i),Form("hDVar_%d_76",i),200,0.,200.);
    hDVar[i][77]= new TH1D(Form("hDVar_%d_77",i),Form("hDVar_%d_77",i),180,0.,180.);
    hDVar[i][78]= new TH1D(Form("hDVar_%d_78",i),Form("hDVar_%d_78",i),180,0.,180.);
    hDVar[i][79]= new TH1D(Form("hDVar_%d_79",i),Form("hDVar_%d_79",i),160,0.,160.);
    hDVar[i][80]= new TH1D(Form("hDVar_%d_80",i),Form("hDVar_%d_80",i),160,0.,160.);

    hDVar[i][81]= new TH1D(Form("hDVar_%d_81",i),Form("hDVar_%d_81",i),200,0.,200.);
    hDVar[i][82]= new TH1D(Form("hDVar_%d_82",i),Form("hDVar_%d_82",i),100,-0.5,99.5);
    hDVar[i][83]= new TH1D(Form("hDVar_%d_83",i),Form("hDVar_%d_83",i),200,0.,200.);
    hDVar[i][84]= new TH1D(Form("hDVar_%d_84",i),Form("hDVar_%d_84",i),200,-10.,10.);
    hDVar[i][85]= new TH1D(Form("hDVar_%d_85",i),Form("hDVar_%d_85",i),200,-10.,10.);
    hDVar[i][86]= new TH1D(Form("hDVar_%d_86",i),Form("hDVar_%d_86",i),180,0.,180.);
    hDVar[i][87]= new TH1D(Form("hDVar_%d_87",i),Form("hDVar_%d_87",i),180,0.,180.);
    hDVar[i][88]= new TH1D(Form("hDVar_%d_88",i),Form("hDVar_%d_88",i),200,0.,200.);
    hDVar[i][89]= new TH1D(Form("hDVar_%d_89",i),Form("hDVar_%d_89",i),160,0.,160.);

    hDVar[i][90]= new TH1D(Form("hDVar_%d_90",i),Form("hDVar_%d_90",i),180,0.,180.);
    hDVar[i][91]= new TH1D(Form("hDVar_%d_91",i),Form("hDVar_%d_91",i),180,0.,180.);
    hDVar[i][92]= new TH1D(Form("hDVar_%d_92",i),Form("hDVar_%d_92",i),180,0.,180.);
    hDVar[i][93]= new TH1D(Form("hDVar_%d_93",i),Form("hDVar_%d_93",i),180,0.,180.);
    hDVar[i][94]= new TH1D(Form("hDVar_%d_94",i),Form("hDVar_%d_94",i),180,0.,180.);
    hDVar[i][95]= new TH1D(Form("hDVar_%d_95",i),Form("hDVar_%d_95",i),180,0.,180.);
  }
  Double_t	  processID;
  Double_t	  mu1IsoCH;
  Double_t	  mu1IsoEM;
  Double_t	  mu1IsoNH;
  Double_t	  mu2IsoCH;
  Double_t	  mu2IsoEM;
  Double_t	  mu2IsoNH;
  Double_t	  mu3IsoCH;
  Double_t	  mu3IsoEM;
  Double_t	  mu3IsoNH;
  Double_t	  mu4IsoCH;
  Double_t	  mu4IsoEM;
  Double_t	  mu4IsoNH;
  Double_t	  ele1IsoCH;
  Double_t	  ele1IsoEM;
  Double_t	  ele1IsoNH;
  Double_t	  ele2IsoCH;
  Double_t	  ele2IsoEM;
  Double_t	  ele2IsoNH;
  Double_t	  ele3IsoCH;
  Double_t	  ele3IsoEM;
  Double_t	  ele3IsoNH;
  Double_t	  ele4IsoCH;
  Double_t	  ele4IsoEM;
  Double_t	  ele4IsoNH;
  Double_t	  chargeMulti;
  Double_t	  visibleE;
  Double_t	  type;
  Double_t	  lq1;
  Double_t	  lq2;
  Double_t	  lq3;
  Double_t	  lq4;
  TLorentzVector*  missingE = 0;
  TLorentzVector*  lepton1 = 0;
  TLorentzVector*  lepton2 = 0;
  TLorentzVector*  lepton3 = 0;
  TLorentzVector*  lepton4 = 0;
  TLorentzVector*  jet1 = 0;
  TLorentzVector*  jet2 = 0;
  TLorentzVector*  jet3 = 0;
  TLorentzVector*  jet4 = 0;
  TLorentzVector*  jet5 = 0;
  TLorentzVector*  jet6 = 0;
  Double_t	   jet1Btag;
  Double_t	   jet2Btag;
  Double_t	   jet3Btag;
  Double_t	   jet4Btag;
  Double_t	   jet5Btag;
  Double_t	   jet6Btag;

  double jetPMin =  5.0;
  double lepPMin = 10.0;
  double ECM      = 240.0;
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
  //for(UInt_t ifile=0; ifile<1; ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    TFile *file = TFile::Open(infilenamev[ifile]);
    file->cd();
    TTree *tree = dynamic_cast<TTree*>(file->Get("lep3Tree"));
    processID = 0;
    tree->SetBranchAddress( "processID"    , &processID   );
    tree->SetBranchAddress( "mu1IsoCH"     , &mu1IsoCH    );
    tree->SetBranchAddress( "mu1IsoEM"     , &mu1IsoEM    );
    tree->SetBranchAddress( "mu1IsoNH"     , &mu1IsoNH    );
    tree->SetBranchAddress( "mu2IsoCH"     , &mu2IsoCH    );
    tree->SetBranchAddress( "mu2IsoEM"     , &mu2IsoEM    );
    tree->SetBranchAddress( "mu2IsoNH"     , &mu2IsoNH    );
    tree->SetBranchAddress( "mu3IsoCH"     , &mu3IsoCH    );
    tree->SetBranchAddress( "mu3IsoEM"     , &mu3IsoEM    );
    tree->SetBranchAddress( "mu3IsoNH"     , &mu3IsoNH    );
    tree->SetBranchAddress( "mu4IsoCH"     , &mu4IsoCH    );
    tree->SetBranchAddress( "mu4IsoEM"     , &mu4IsoEM    );
    tree->SetBranchAddress( "mu4IsoNH"     , &mu4IsoNH    );
    tree->SetBranchAddress( "ele1IsoCH"    , &ele1IsoCH   );
    tree->SetBranchAddress( "ele1IsoEM"    , &ele1IsoEM   );
    tree->SetBranchAddress( "ele1IsoNH"    , &ele1IsoNH   );
    tree->SetBranchAddress( "ele2IsoCH"    , &ele2IsoCH   );
    tree->SetBranchAddress( "ele2IsoEM"    , &ele2IsoEM   );
    tree->SetBranchAddress( "ele2IsoNH"    , &ele2IsoNH   );
    tree->SetBranchAddress( "ele3IsoCH"    , &ele3IsoCH   );
    tree->SetBranchAddress( "ele3IsoEM"    , &ele3IsoEM   );
    tree->SetBranchAddress( "ele3IsoNH"    , &ele3IsoNH   );
    tree->SetBranchAddress( "ele4IsoCH"    , &ele4IsoCH   );
    tree->SetBranchAddress( "ele4IsoEM"    , &ele4IsoEM   );
    tree->SetBranchAddress( "ele4IsoNH"    , &ele4IsoNH   );
    tree->SetBranchAddress( "chargeMulti"  , &chargeMulti );
    tree->SetBranchAddress( "visibleE"	   , &visibleE    );
    tree->SetBranchAddress( "type"	   , &type	  );
    tree->SetBranchAddress( "lq1"	   , &lq1	  );
    tree->SetBranchAddress( "lq2"	   , &lq2	  );
    tree->SetBranchAddress( "lq3"	   , &lq3	  );
    tree->SetBranchAddress( "lq4"	   , &lq4	  );
    tree->SetBranchAddress( "missingE"     , &missingE    );
    tree->SetBranchAddress( "lepton1"	   , &lepton1	  );
    tree->SetBranchAddress( "lepton2"	   , &lepton2	  );
    tree->SetBranchAddress( "lepton3"	   , &lepton3	  );
    tree->SetBranchAddress( "lepton4"	   , &lepton4	  );
    tree->SetBranchAddress( "jet1"	   , &jet1	  );
    tree->SetBranchAddress( "jet2"	   , &jet2	  );
    tree->SetBranchAddress( "jet3"	   , &jet3	  );
    tree->SetBranchAddress( "jet4"	   , &jet4	  );
    tree->SetBranchAddress( "jet5"	   , &jet5	  );
    tree->SetBranchAddress( "jet6"	   , &jet6	  );
    tree->SetBranchAddress( "jet1Btag"     , &jet1Btag    );
    tree->SetBranchAddress( "jet2Btag"     , &jet2Btag    );
    tree->SetBranchAddress( "jet3Btag"     , &jet3Btag    );
    tree->SetBranchAddress( "jet4Btag"     , &jet4Btag    );
    tree->SetBranchAddress( "jet5Btag"     , &jet5Btag    );
    tree->SetBranchAddress( "jet6Btag"     , &jet6Btag    );

    TH1D* hDNormalization = (TH1D*)(file->Get("hNormalization"));
    Double_t nTotEvt = hDNormalization->GetSumOfWeights();
    Double_t weight  = scale1fb[ifile] * lumi / nTotEvt;
    UInt_t NEvents = tree->GetEntries(); //if(ifile == 0) {NEvents = NEvents/5; weight = weight*5.0;};
    printf("cross-section: %f lumi: %f nTotEvt: %f --> weight: %f | NEvents: %d\n",scale1fb[ifile],lumi,nTotEvt,weight,NEvents);
    for(UInt_t ientry = 0; ientry <NEvents; ientry++){
      tree->GetEntry(ientry);
      if (ientry % 100000 == 0) printf("Event %d of %d\n",ientry,(int)NEvents);

      Int_t iHist  = ifile;
      Int_t iClass = 0;
      if     (ifile == 0 && processID != 24) {
        if     (processID ==  5) {iHist = infilenamev.size()+0;iClass = 2;}
        else if(processID ==  4) {iHist = infilenamev.size()+1;iClass = 3;}
        else if(processID ==  3) {iHist = infilenamev.size()+2;iClass = 4;}
        else if(processID == 15) {iHist = infilenamev.size()+3;iClass = 5;}
        else if(processID == 13) {iHist = infilenamev.size()+4;iClass = 6;}
        else if(processID == 21) {iHist = infilenamev.size()+5;iClass = 7;}
        else if(processID == 22) {iHist = infilenamev.size()+6;iClass = 8;}
        else if(processID == 23) {iHist = infilenamev.size()+7;iClass = 9;}
	else {cout << "forbidden decay: " << processID << endl; assert(0);}
      }
      else if(ifile == 0 && processID == 24) {
        iClass = 1;
      }

      hDVar[iHist][0]->Fill(TMath::Min(chargeMulti,99.499),weight);
      if(chargeMulti <= 1) continue;
      hDVar[iHist][1]->Fill(TMath::Min(visibleE,299.999),weight);
      if(visibleE <= 50) continue;

      if(ifile == 0) hDSig[0]->Fill(processID,weight);

      int totalQ = 0;
      // Lepton selection
      int lType[4] = {(int)type%10*lq1,(int)(((int)type%100)/10)*lq2,(int)(type/100)%10*lq3,(int)type/1000*lq4};
      TLorentzVector* leptonG[4] = {0,0,0,0};
      if(lepton1->P() > lepPMin && ((abs(lType[0]) == 1 && (mu1IsoCH+mu1IsoEM+mu1IsoNH   ) < 5.0) ||
                   	            (abs(lType[0]) == 2 && (ele1IsoCH+ele1IsoEM+ele1IsoNH) < 5.0))) leptonG[0] = lepton1;
      if(lepton2->P() > lepPMin && ((abs(lType[1]) == 1 && (mu2IsoCH+mu2IsoEM+mu2IsoNH   ) < 5.0) ||
                   	            (abs(lType[1]) == 2 && (ele2IsoCH+ele2IsoEM+ele2IsoNH) < 5.0))) leptonG[1] = lepton2;
      if(lepton3->P() > lepPMin && ((abs(lType[2]) == 1 && (mu3IsoCH+mu3IsoEM+mu3IsoNH   ) < 5.0) ||
                   	            (abs(lType[2]) == 2 && (ele3IsoCH+ele3IsoEM+ele3IsoNH) < 5.0))) leptonG[2] = lepton3;
      if(lepton4->P() > lepPMin && ((abs(lType[3]) == 1 && (mu4IsoCH+mu4IsoEM+mu4IsoNH   ) < 5.0) ||
                            	    (abs(lType[3]) == 2 && (ele4IsoCH+ele4IsoEM+ele4IsoNH) < 5.0))) leptonG[3] = lepton4;
      for(int i=0; i<4; i++) if(leptonG[i] == 0) lType[i] = 0;
      if(leptonG[0] != 0) totalQ += lq1; if(leptonG[1] != 0) totalQ += lq2; 
      if(leptonG[2] != 0) totalQ += lq3; if(leptonG[3] != 0) totalQ += lq4; 
      totalQ = TMath::Abs(totalQ);

      // Jet selection
      Bool_t badJet[6] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(badJet[0] == kFALSE && ((leptonG[0] !=  0 && jet1->P() > 0 && leptonG[0]->DeltaR(*jet1) < 0.3)||jet1->P() < jetPMin)) badJet[0] = kTRUE;
      if(badJet[1] == kFALSE && ((leptonG[0] !=  0 && jet2->P() > 0 && leptonG[0]->DeltaR(*jet2) < 0.3)||jet2->P() < jetPMin)) badJet[1] = kTRUE;
      if(badJet[2] == kFALSE && ((leptonG[0] !=  0 && jet3->P() > 0 && leptonG[0]->DeltaR(*jet3) < 0.3)||jet3->P() < jetPMin)) badJet[2] = kTRUE;
      if(badJet[3] == kFALSE && ((leptonG[0] !=  0 && jet4->P() > 0 && leptonG[0]->DeltaR(*jet4) < 0.3)||jet4->P() < jetPMin)) badJet[3] = kTRUE;
      if(badJet[4] == kFALSE && ((leptonG[0] !=  0 && jet5->P() > 0 && leptonG[0]->DeltaR(*jet5) < 0.3)||jet5->P() < jetPMin)) badJet[4] = kTRUE;
      if(badJet[5] == kFALSE && ((leptonG[0] !=  0 && jet6->P() > 0 && leptonG[0]->DeltaR(*jet6) < 0.3)||jet6->P() < jetPMin)) badJet[5] = kTRUE;
      if(badJet[0] == kFALSE && ((leptonG[1] !=  0 && jet1->P() > 0 && leptonG[1]->DeltaR(*jet1) < 0.3)||jet1->P() < jetPMin)) badJet[0] = kTRUE;
      if(badJet[1] == kFALSE && ((leptonG[1] !=  0 && jet2->P() > 0 && leptonG[1]->DeltaR(*jet2) < 0.3)||jet2->P() < jetPMin)) badJet[1] = kTRUE;
      if(badJet[2] == kFALSE && ((leptonG[1] !=  0 && jet3->P() > 0 && leptonG[1]->DeltaR(*jet3) < 0.3)||jet3->P() < jetPMin)) badJet[2] = kTRUE;
      if(badJet[3] == kFALSE && ((leptonG[1] !=  0 && jet4->P() > 0 && leptonG[1]->DeltaR(*jet4) < 0.3)||jet4->P() < jetPMin)) badJet[3] = kTRUE;
      if(badJet[4] == kFALSE && ((leptonG[1] !=  0 && jet5->P() > 0 && leptonG[1]->DeltaR(*jet5) < 0.3)||jet5->P() < jetPMin)) badJet[4] = kTRUE;
      if(badJet[5] == kFALSE && ((leptonG[1] !=  0 && jet6->P() > 0 && leptonG[1]->DeltaR(*jet6) < 0.3)||jet6->P() < jetPMin)) badJet[5] = kTRUE;
      if(badJet[0] == kFALSE && ((leptonG[2] !=  0 && jet1->P() > 0 && leptonG[2]->DeltaR(*jet1) < 0.3)||jet1->P() < jetPMin)) badJet[0] = kTRUE;
      if(badJet[1] == kFALSE && ((leptonG[2] !=  0 && jet2->P() > 0 && leptonG[2]->DeltaR(*jet2) < 0.3)||jet2->P() < jetPMin)) badJet[1] = kTRUE;
      if(badJet[2] == kFALSE && ((leptonG[2] !=  0 && jet3->P() > 0 && leptonG[2]->DeltaR(*jet3) < 0.3)||jet3->P() < jetPMin)) badJet[2] = kTRUE;
      if(badJet[3] == kFALSE && ((leptonG[2] !=  0 && jet4->P() > 0 && leptonG[2]->DeltaR(*jet4) < 0.3)||jet4->P() < jetPMin)) badJet[3] = kTRUE;
      if(badJet[4] == kFALSE && ((leptonG[2] !=  0 && jet5->P() > 0 && leptonG[2]->DeltaR(*jet5) < 0.3)||jet5->P() < jetPMin)) badJet[4] = kTRUE;
      if(badJet[5] == kFALSE && ((leptonG[2] !=  0 && jet6->P() > 0 && leptonG[2]->DeltaR(*jet6) < 0.3)||jet6->P() < jetPMin)) badJet[5] = kTRUE;
      if(badJet[0] == kFALSE && ((leptonG[3] !=  0 && jet1->P() > 0 && leptonG[3]->DeltaR(*jet1) < 0.3)||jet1->P() < jetPMin)) badJet[0] = kTRUE;
      if(badJet[1] == kFALSE && ((leptonG[3] !=  0 && jet2->P() > 0 && leptonG[3]->DeltaR(*jet2) < 0.3)||jet2->P() < jetPMin)) badJet[1] = kTRUE;
      if(badJet[2] == kFALSE && ((leptonG[3] !=  0 && jet3->P() > 0 && leptonG[3]->DeltaR(*jet3) < 0.3)||jet3->P() < jetPMin)) badJet[2] = kTRUE;
      if(badJet[3] == kFALSE && ((leptonG[3] !=  0 && jet4->P() > 0 && leptonG[3]->DeltaR(*jet4) < 0.3)||jet4->P() < jetPMin)) badJet[3] = kTRUE;
      if(badJet[4] == kFALSE && ((leptonG[3] !=  0 && jet5->P() > 0 && leptonG[3]->DeltaR(*jet5) < 0.3)||jet5->P() < jetPMin)) badJet[4] = kTRUE;
      if(badJet[5] == kFALSE && ((leptonG[3] !=  0 && jet6->P() > 0 && leptonG[3]->DeltaR(*jet6) < 0.3)||jet6->P() < jetPMin)) badJet[5] = kTRUE;
      TLorentzVector*  jetG[6] = {0,0,0,0,0,0};
      if     (badJet[0] == kFALSE                ) jetG[0] = jet1;
      if     (badJet[1] == kFALSE && jetG[0] == 0) jetG[0] = jet2;
      else if(badJet[1] == kFALSE                ) jetG[1] = jet2;
      if     (badJet[2] == kFALSE && jetG[0] == 0) jetG[0] = jet3;
      else if(badJet[2] == kFALSE && jetG[1] == 0) jetG[1] = jet3;
      else if(badJet[2] == kFALSE                ) jetG[2] = jet3;
      if     (badJet[3] == kFALSE && jetG[0] == 0) jetG[0] = jet4;
      else if(badJet[3] == kFALSE && jetG[1] == 0) jetG[1] = jet4;
      else if(badJet[3] == kFALSE && jetG[2] == 0) jetG[2] = jet4;
      else if(badJet[3] == kFALSE                ) jetG[3] = jet4;
      if     (badJet[4] == kFALSE && jetG[0] == 0) jetG[0] = jet5;
      else if(badJet[4] == kFALSE && jetG[1] == 0) jetG[1] = jet5;
      else if(badJet[4] == kFALSE && jetG[2] == 0) jetG[2] = jet5;
      else if(badJet[4] == kFALSE && jetG[3] == 0) jetG[3] = jet5;
      else if(badJet[4] == kFALSE                ) jetG[4] = jet5;
      if     (badJet[5] == kFALSE && jetG[0] == 0) jetG[0] = jet6;
      else if(badJet[5] == kFALSE && jetG[1] == 0) jetG[1] = jet6;
      else if(badJet[5] == kFALSE && jetG[2] == 0) jetG[2] = jet6;
      else if(badJet[5] == kFALSE && jetG[3] == 0) jetG[3] = jet6;
      else if(badJet[5] == kFALSE && jetG[4] == 0) jetG[4] = jet6;
      else if(badJet[5] == kFALSE                ) jetG[5] = jet6;
      
      double mTotV[4] = {0,0,0,0}; double mTot = 0.0;
      int nLeptons = 0; for(int i=0; i<4; i++) if(leptonG[i] != 0) {nLeptons++;mTotV[0]+=leptonG[i]->Px();mTotV[1]+=leptonG[i]->Py();mTotV[2]+=leptonG[i]->Pz();mTotV[3]+=leptonG[i]->P();}
      int nJets    = 0; for(int i=0; i<6; i++) if(jetG[i] != 0) {nJets++;mTotV[0]+=jetG[i]->Px();mTotV[1]+=jetG[i]->Py();mTotV[2]+=jetG[i]->Pz();mTotV[3]+=jetG[i]->P();}
      mTot = mTotV[3]*mTotV[3]-mTotV[0]*mTotV[0]-mTotV[1]*mTotV[1]-mTotV[2]*mTotV[2];
      if(mTot>0) mTot=sqrt(mTot); else mTot = 0.0;
      hDVar[iHist][2]->Fill((double)(nJets+10*nLeptons),weight);
      if(nLeptons == 2) hDVar[iHist][3]->Fill((double)totalQ,weight);
      if(nLeptons == 3) hDVar[iHist][4]->Fill((double)totalQ,weight);
      if(nLeptons == 4) hDVar[iHist][5]->Fill((double)totalQ,weight);

      // Missing energy
      double MET[3] = {0, 0, 0};
      for(int i=0; i<4; i++) if(leptonG[i] != 0) {MET[0] -= leptonG[i]->Px();MET[1] -= leptonG[i]->Py();MET[2] -= leptonG[i]->Pz();}
      for(int i=0; i<6; i++) if(jetG[i] != 0) {MET[0] -= jetG[i]->Px();MET[1] -= jetG[i]->Py();MET[2] -= jetG[i]->Pz();}
      double MissingMass = (ECM-visibleE)*(ECM-visibleE)-MET[0]*MET[0]-MET[1]*MET[1]-MET[2]*MET[2];
      if(MissingMass > 0) MissingMass = sqrt(MissingMass); else MissingMass = 0.0;
      hDVar[iHist][6]->Fill(TMath::Min(MissingMass,299.99),weight);

      const TVector3 theMET(MET[0],MET[1],MET[2]);

      // Finding Z->ll decays
      if(verbose) cout << lType[0] << " "  << lType[1] << " "  << lType[2] << " "  << lType[3] << " "  <<type<<endl;
      int iZllCand[2] = {-1,-1};
      double massCloseZll = 1000.;
      for(int i=0; i<4; i++){
        for(int j=i+1; j<4; j++){
	  if(lType[i] != 0 && lType[j] != 0 && lType[i]*lType[j] < 0 && abs(lType[i]) == abs(lType[j])) {
	    if(TMath::Abs((*leptonG[i]+*leptonG[j]).M()-91.1876) < TMath::Abs(massCloseZll-91.1876)){
	      massCloseZll =(*leptonG[i]+*leptonG[j]).M();
	      iZllCand[0] = i;
	      iZllCand[1] = j;
	    }
	  }
	}
      }
      hDVar[iHist][7]->Fill(TMath::Min(massCloseZll,149.99),weight);
      if(TMath::Abs(massCloseZll-91.1876) >= 15) {iZllCand[0] = -1; iZllCand[1] = -1;}

      // Finding Z->qq decays
      int iZqqCand[2] = {-1,-1};
      double massCloseZqq = 1000.;
      for(int i=0; i<6; i++){
        for(int j=i+1; j<6; j++){
	  if(jetG[i] != 0 && jetG[j] != 0) {
	    if(TMath::Abs((*jetG[i]+*jetG[j]).M()-91.1876) < TMath::Abs(massCloseZqq-91.1876)){
	      massCloseZqq =(*jetG[i]+*jetG[j]).M();
	      iZqqCand[0] = i;
	      iZqqCand[1] = j;
	    }
	  }
	}
      }
      hDVar[iHist][8]->Fill(TMath::Min(massCloseZqq,149.99),weight);
      double mRecoil[3] = {0,0,0};
      if(TMath::Abs(massCloseZll-91.1876) < 15){
        double pZ = (*leptonG[iZllCand[0]]+*leptonG[iZllCand[1]]).P();
	double eZ = sqrt(massCloseZll*massCloseZll+pZ*pZ);
	mRecoil[0] = (240.0-eZ)*(240.0-eZ)-pZ*pZ;
	if(mRecoil[0] > 0) mRecoil[0] = sqrt(mRecoil[0]); else mRecoil[0] = 0.0;
	hDVar[iHist][9]->Fill(TMath::Min(mRecoil[0],199.99),weight);
      }
      if(TMath::Abs(massCloseZqq-91.1876) < 25){
        double pZ = (*jetG[iZqqCand[0]]+*jetG[iZqqCand[1]]).P();
	double eZ = sqrt(91.1876*91.1876+pZ*pZ);
	mRecoil[1] = (240.0-eZ)*(240.0-eZ)-pZ*pZ;
	if(mRecoil[1] > 0) mRecoil[1] = sqrt(mRecoil[1]); else mRecoil[1] = 0.0;
	hDVar[iHist][10]->Fill(TMath::Min(mRecoil[1],199.99),weight);
      }
      if(MissingMass > 60.0){
        double pZ = sqrt(missingE->Pt()*missingE->Pt()+theMET.Pz()*theMET.Pz());
	double eZ = sqrt(91.1876*91.1876+pZ*pZ);
	mRecoil[2] = (235.0-eZ)*(235.0-eZ)-pZ*pZ;
	if(mRecoil[2] > 0) mRecoil[2] = sqrt(mRecoil[2]); else mRecoil[2] = 0.0;
	hDVar[iHist][11]->Fill(TMath::Min(mRecoil[2],199.99),weight);
      }

      // Finding W->qq decays
      double massCloseWqq = 1000.;
      double massCloseWqqNoZRemoval = 1000.;
      for(int i=0; i<6; i++){
        for(int j=i+1; j<6; j++){
	  if(jetG[i] != 0 && jetG[j] != 0 && iZqqCand[0] != i && iZqqCand[1] != i 
	                                  && iZqqCand[0] != j && iZqqCand[1] != j) {
	    if(TMath::Abs((*jetG[i]+*jetG[j]).M()-80.40) < TMath::Abs(massCloseWqq-80.40)){
	      massCloseWqq =(*jetG[i]+*jetG[j]).M();
	    }
	  }
	  if(jetG[i] != 0 && jetG[j] != 0) {
	    if(TMath::Abs((*jetG[i]+*jetG[j]).M()-80.40) < TMath::Abs(massCloseWqqNoZRemoval-80.40)){
	      massCloseWqqNoZRemoval =(*jetG[i]+*jetG[j]).M();
	    }
	  }
	}
      }
      hDVar[iHist][12]->Fill(TMath::Min(massCloseWqqNoZRemoval,149.99),weight);

      double jetBtagMax[3] = {-1000.0,-1000.0,-1000.0};
      if     (jet1Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet1Btag;}
      else if(jet1Btag > jetBtagMax[1]) {jetBtagMax[1] = jet1Btag;}
      if     (jet2Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet2Btag;}
      else if(jet2Btag > jetBtagMax[1]) {jetBtagMax[1] = jet2Btag;}
      if     (jet3Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet3Btag;}
      else if(jet3Btag > jetBtagMax[1]) {jetBtagMax[1] = jet3Btag;}
      if     (jet4Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet4Btag;}
      else if(jet4Btag > jetBtagMax[1]) {jetBtagMax[1] = jet4Btag;}
      if     (jet5Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet5Btag;}
      else if(jet5Btag > jetBtagMax[1]) {jetBtagMax[1] = jet5Btag;}
      if     (jet6Btag > jetBtagMax[0]) {jetBtagMax[1] = jetBtagMax[0]; jetBtagMax[0] = jet6Btag;}
      else if(jet6Btag > jetBtagMax[1]) {jetBtagMax[1] = jet6Btag;}

      if     (jet1Btag > jetBtagMax[2] && iZqqCand[0] != 1 && iZqqCand[1] != 1) {jetBtagMax[2] = jet1Btag;}
      if     (jet2Btag > jetBtagMax[2] && iZqqCand[0] != 2 && iZqqCand[1] != 2) {jetBtagMax[2] = jet2Btag;}
      if     (jet3Btag > jetBtagMax[2] && iZqqCand[0] != 3 && iZqqCand[1] != 3) {jetBtagMax[2] = jet3Btag;}
      if     (jet4Btag > jetBtagMax[2] && iZqqCand[0] != 4 && iZqqCand[1] != 4) {jetBtagMax[2] = jet4Btag;}
      if     (jet5Btag > jetBtagMax[2] && iZqqCand[0] != 5 && iZqqCand[1] != 5) {jetBtagMax[2] = jet5Btag;}
      if     (jet6Btag > jetBtagMax[2] && iZqqCand[0] != 6 && iZqqCand[1] != 6) {jetBtagMax[2] = jet6Btag;}

      for(int i=0; i<3; i++) jetBtagMax[i] = TMath::Max(TMath::Min(jetBtagMax[i],9.9999),-9.9999);

      double angleLMET = 1000;
      double anglesQQ[2] = {0.,1000.}; double angleQL = 1000; double etaJetMax = 0;
      for(int i=0; i<6; i++){
        if(jetG[i] != 0) if(TMath::Abs(jetG[i]->Eta()) > etaJetMax) etaJetMax = TMath::Abs(jetG[i]->Eta());
        for(int j=i+1; j<6; j++) {
          if(jetG[i] != 0 && jetG[j] != 0){
            const TVector3 a(jetG[j]->Px(),jetG[j]->Py(),jetG[j]->Pz());
            double angle = jetG[i]->Angle(a);
	    if(angle > anglesQQ[0]) anglesQQ[0] = angle;
	    if(angle < anglesQQ[1]) anglesQQ[1] = angle;
          }
        }
        for(int j=0; j<4; j++){
          if(jetG[i] != 0 && leptonG[j] != 0){
            const TVector3 a(leptonG[j]->Px(),leptonG[j]->Py(),leptonG[j]->Pz());
            double angle = jetG[i]->Angle(a);
	    if(angle < angleQL) angleQL = angle;
          }
          if(leptonG[j] != 0){
            const TVector3 a(leptonG[j]->Px(),leptonG[j]->Py(),leptonG[j]->Pz());
	    double angle = theMET.Angle(a);
	    if(angle < angleLMET) angleLMET = angle;
	  }
	}
      }
      anglesQQ[0] = anglesQQ[0]*180/TMath::Pi();
      anglesQQ[1] = anglesQQ[1]*180/TMath::Pi();
      angleQL     = angleQL    *180/TMath::Pi();

      // Projected MET
      double PMET = theMET.Mag();
      if(angleLMET < TMath::Pi()/2) PMET = PMET * sin(angleLMET);
      angleLMET   = angleLMET  *180/TMath::Pi();
      
      // Final selection

      // ZH->qqqqqq
      if(MissingMass < 60 && nJets >= 6 && nLeptons == 0 && mRecoil[1] > 0){
        hDVar[iHist][13]->Fill(TMath::Min(etaJetMax,2.999),weight);
        hDVar[iHist][14]->Fill(TMath::Min(chargeMulti,99.499),weight);
        hDVar[iHist][15]->Fill(TMath::Min(mTot,299.999),weight);
	if(chargeMulti > 35 && mTot > 140 && etaJetMax < 2.5 && jetBtagMax[2] < 4.0){	  
	  int indexW[4] = {-1,-1,-1,-1}; int theIndex = 0;
	  for(int i=0; i<6; i++) if(jetG[i] != 0 && iZqqCand[0] != i && iZqqCand[1] != i) {
	    indexW[theIndex] = i;
	    theIndex++;
	  }
          double massW1 = (*jetG[indexW[0]]+*jetG[indexW[1]]).M();
	  double massW2 = (*jetG[indexW[2]]+*jetG[indexW[3]]).M();
          hDVar[iHist][16]->Fill(TMath::Min(massCloseZqq,149.99),weight);
	  hDVar[iHist][17]->Fill(TMath::Min(mRecoil[1],159.999),weight);
	  hDVar[iHist][18]->Fill(TMath::Min(massW1,199.99),weight);
	  hDVar[iHist][19]->Fill(TMath::Min(massW2,199.99),weight);	  
	  hDVar[iHist][20]->Fill(TMath::Min(theMET.Pt(),199.99),weight);
          hDVar[iHist][21]->Fill(TMath::Min(TMath::Abs(theMET.Pz()),199.999),weight);
          hDVar[iHist][22]->Fill(jetBtagMax[0],weight);
          hDVar[iHist][23]->Fill(jetBtagMax[2],weight);
	  if(ifile == 0) hDSig[1]->Fill(processID,weight);
	  if(mRecoil[1] > 120){
	    nClass[0][iClass] += weight;
	  }
	}
      }

      // ZH->qqqqln
      if(nLeptons == 1 && nJets >= 4 && mRecoil[1] > 0 && jetBtagMax[2] < 4.0){
	double pLep = 0.0;
	for(int i=0; i<4; i++) if(leptonG[i] != 0) pLep = leptonG[i]->E();
	hDVar[iHist][24]->Fill(TMath::Min(pLep,199.99),weight);
        hDVar[iHist][25]->Fill(anglesQQ[0],weight);
        hDVar[iHist][26]->Fill(TMath::Min(PMET,179.99),weight);
	if(chargeMulti > 20 && pLep < 55.0 && anglesQQ[0] > 130 && 
	   PMET > 25.0){
          hDVar[iHist][27]->Fill(TMath::Min(massCloseZqq,199.99),weight);
	  hDVar[iHist][28]->Fill(TMath::Min(theMET.Pt(),199.99),weight);
          hDVar[iHist][29]->Fill(TMath::Min(massCloseWqq,199.999),weight);
          hDVar[iHist][30]->Fill(TMath::Min(mRecoil[1],159.99),weight);
	  if(ifile == 0) hDSig[2]->Fill(processID,weight);
          hDVar[iHist][90]->Fill(angleLMET,weight);
          hDVar[iHist][96]->Fill(jetBtagMax[0],weight);
          hDVar[iHist][97]->Fill(jetBtagMax[2],weight);
	  
	  if(mRecoil[1] > 120){
	    nClass[1][iClass] += weight;
	  }
	}
      }

      int indexW[2];int index = 0;
      for(int i=0; i<4; i++) if(leptonG[i] != 0 && 
          iZllCand[0] != i && iZllCand[1] != i) {indexW[index] = i; index++;}
      // ZH->qqlnln
      if(nLeptons == 2 && totalQ == 0 && nJets == 2 && iZllCand[0] == -1 && mRecoil[1] > 0){
	double phiLL = DeltaPhi(leptonG[indexW[0]]->Phi(),leptonG[indexW[1]]->Phi())*180/TMath::Pi();
        hDVar[iHist][31]->Fill(TMath::Min(chargeMulti,99.499),weight);
	hDVar[iHist][32]->Fill(TMath::Min(TMath::Max(leptonG[indexW[0]]->E(),leptonG[indexW[1]]->E()),99.99),weight);	
	hDVar[iHist][33]->Fill(TMath::Min(phiLL,199.99),weight);
	hDVar[iHist][34]->Fill(TMath::Min(PMET,199.99),weight);
        hDVar[iHist][35]->Fill(TMath::Min((*leptonG[indexW[0]]+*leptonG[indexW[1]]).M(),199.999),weight);
 	if(chargeMulti > 10 && (*leptonG[indexW[0]]+*leptonG[indexW[1]]).M() < 70 &&
	                       (*leptonG[indexW[0]]+*leptonG[indexW[1]]).M() > 10 &&
           TMath::Max(leptonG[indexW[0]]->E(),leptonG[indexW[1]]->E()) < 70 &&
	   PMET > 25){
	  if(ifile == 0) hDSig[3]->Fill(processID,weight);
          hDVar[iHist][91]->Fill(angleLMET,weight);
          if(abs(lType[indexW[0]]) != abs(lType[indexW[1]])) {
	    hDVar[iHist][36]->Fill(TMath::Min(mRecoil[1],159.99),weight);
	    if(mRecoil[1] > 120){
	      nClass[2][iClass] += weight;
	    }
	  } 
	  else {
	    hDVar[iHist][37]->Fill(TMath::Min(mRecoil[1],159.99),weight);
	    if(mRecoil[1] > 120){
	      nClass[3][iClass] += weight;
	    }
	  }
        }
      }

      // ZH->nnqqqq
      if(MissingMass > 40 && nJets >= 4 && nJets <= 5 && nLeptons == 0 && mRecoil[2] > 0 && jetBtagMax[0] < 4.0){
        hDVar[iHist][38]->Fill(TMath::Min(MissingMass,199.99),weight);
        hDVar[iHist][39]->Fill(TMath::Min(chargeMulti,99.499),weight);
	hDVar[iHist][40]->Fill(TMath::Min(theMET.Pt(),199.99),weight);
        hDVar[iHist][41]->Fill(TMath::Min(TMath::Abs(theMET.Pz()),99.999),weight);
	if(MissingMass > 70 && MissingMass < 130 && chargeMulti > 15 && TMath::Abs(theMET.Pz()) < 30.0 && 
	   theMET.Pt() > 20.0 && etaJetMax < 2.5){
	  int indexWqq[4] = {-1,-1,-1,-1}; int theIndex = 0;
	  for(int i=0; i<6; i++) if(jetG[i] != 0) {
	    indexWqq[theIndex] = i;
	    theIndex++;
	  }
          double massW1 = (*jetG[indexWqq[0]]+*jetG[indexWqq[1]]).M();
	  double massW2 = (*jetG[indexWqq[2]]+*jetG[indexWqq[3]]).M();
	  hDVar[iHist][42]->Fill(TMath::Min(massW1,199.99),weight);
	  hDVar[iHist][43]->Fill(TMath::Min(massW2,199.99),weight);	  
          hDVar[iHist][44]->Fill(jetBtagMax[0],weight);
          hDVar[iHist][45]->Fill(jetBtagMax[1],weight);
	  hDVar[iHist][46]->Fill(TMath::Min(mRecoil[2],159.99),weight);	  
	  if(ifile == 0) hDSig[4]->Fill(processID,weight);
	  if(mRecoil[2] > 115){
	    nClass[4][iClass] += weight;
	  }
        }
      }

      // ZH->nnlnln
      if(nLeptons == 2 && totalQ == 0 && nJets == 0 && iZllCand[0] == -1 && chargeMulti == 2){
	double phiLL = DeltaPhi(leptonG[indexW[0]]->Phi(),leptonG[indexW[1]]->Phi())*180/TMath::Pi();
	hDVar[iHist][47]->Fill(phiLL,weight);
	hDVar[iHist][48]->Fill(TMath::Min(PMET,199.99),weight);
        hDVar[iHist][49]->Fill(TMath::Min(TMath::Max(leptonG[indexW[0]]->E(),leptonG[indexW[1]]->E()),199.999),weight);
        hDVar[iHist][50]->Fill(TMath::Min(TMath::Min(leptonG[indexW[0]]->E(),leptonG[indexW[1]]->E()),199.999),weight);
	hDVar[iHist][51]->Fill(TMath::Min(MissingMass,199.99),weight);
	hDVar[iHist][52]->Fill(TMath::ACos(theMET.CosTheta())*180/TMath::Pi(),weight);
        hDVar[iHist][53]->Fill(TMath::Min((*leptonG[indexW[0]]+*leptonG[indexW[1]]).M(),199.999),weight);
	if(MissingMass > 145 && (*leptonG[indexW[0]]+*leptonG[indexW[1]]).M() < 80 && 
	                        (*leptonG[indexW[0]]+*leptonG[indexW[1]]).M() > 10 && 
	   PMET > 20 && phiLL < 165){
	  if(ifile == 0) hDSig[5]->Fill(processID,weight);
          hDVar[iHist][92]->Fill(angleLMET,weight);
          if(abs(lType[indexW[0]]) != abs(lType[indexW[1]])) {
            hDVar[iHist][54]->Fill(TMath::Min(TMath::Abs(theMET.Pz()),199.999),weight);
	    nClass[5][iClass] += weight;
          }
	  else {
	    hDVar[iHist][55]->Fill(TMath::Min(TMath::Abs(theMET.Pz()),199.999),weight);
	    nClass[6][iClass] += weight;
	  }
	}
      }

      // ZH->llqqqq
      if(nLeptons == 2 && totalQ == 0 && nJets >= 4 &&
         iZllCand[0] != -1 && mRecoil[0] > 0 && jetBtagMax[0] < 4.0){
        hDVar[iHist][56]->Fill(TMath::Min(massCloseZll,199.99),weight);
        hDVar[iHist][57]->Fill(TMath::Min(chargeMulti,99.499),weight);
	hDVar[iHist][58]->Fill(TMath::Min(theMET.Pt(),199.99),weight);
        hDVar[iHist][59]->Fill(TMath::Min(leptonG[iZllCand[0]]->E()+leptonG[iZllCand[1]]->E(),199.999),weight);
        const TVector3 a(leptonG[iZllCand[0]]->Px(),leptonG[iZllCand[0]]->Py(),leptonG[iZllCand[0]]->Pz());
        double angleLL = leptonG[iZllCand[1]]->Angle(a)*180/TMath::Pi();
	hDVar[iHist][60]->Fill(angleLL,weight); 	
	if(chargeMulti > 10 && angleLL > 120 && etaJetMax < 2.5){
          hDVar[iHist][61]->Fill(jetBtagMax[0],weight);
          hDVar[iHist][62]->Fill(jetBtagMax[1],weight);
	  hDVar[iHist][63]->Fill(TMath::Min(mRecoil[0],159.99),weight);	
	  if(ifile == 0) hDSig[6]->Fill(processID,weight);
	  if(mRecoil[0] > 120 && mRecoil[0] < 140){
	    nClass[7][iClass] += weight;
	  }
        }
      }

      // ZH->lllnqq
      if(nLeptons == 3 && totalQ == 1 && nJets >= 1 && chargeMulti > 5 && 
         iZllCand[0] != -1 && mRecoil[0] > 0){
        hDVar[iHist][64]->Fill(TMath::Min(massCloseZll,199.99),weight);
        hDVar[iHist][65]->Fill(TMath::Min(leptonG[iZllCand[0]]->E()+leptonG[iZllCand[1]]->E(),199.999),weight);
        const TVector3 a(leptonG[iZllCand[0]]->Px(),leptonG[iZllCand[0]]->Py(),leptonG[iZllCand[0]]->Pz());
        double angleLL = leptonG[iZllCand[1]]->Angle(a)*180/TMath::Pi();
	hDVar[iHist][66]->Fill(angleLL,weight);
        hDVar[iHist][67]->Fill(TMath::Min(chargeMulti,99.499),weight);
    	hDVar[iHist][68]->Fill(TMath::Min(leptonG[indexW[0]]->P(),199.99),weight);
	hDVar[iHist][69]->Fill(TMath::Min(PMET,199.99),weight);
        hDVar[iHist][98]->Fill(jetBtagMax[0],weight);
        hDVar[iHist][99]->Fill(jetBtagMax[1],weight);
	if(chargeMulti > 10 && angleLL > 120 && leptonG[indexW[0]]->P() < 50 && PMET > 20){
	  hDVar[iHist][70]->Fill(TMath::Min(massCloseWqqNoZRemoval,199.99),weight);
          hDVar[iHist][71]->Fill(angleQL,weight);
          hDVar[iHist][72]->Fill(TMath::Min(mRecoil[0],159.999),weight);
	  if(ifile == 0) hDSig[7]->Fill(processID,weight);
          hDVar[iHist][93]->Fill(angleLMET,weight);
	  if(mRecoil[0] > 120 && mRecoil[0] < 140){
	    nClass[8][iClass] += weight;
	  }
	}
      }

      // ZH->lllnln
      if(nLeptons == 4 && totalQ == 0 && nJets == 0 && 
         iZllCand[0] != -1 && mRecoil[0] > 0){
	double phiLL = DeltaPhi(leptonG[indexW[0]]->Phi(),leptonG[indexW[1]]->Phi())*180/TMath::Pi();
        const TVector3 a(leptonG[iZllCand[0]]->Px(),leptonG[iZllCand[0]]->Py(),leptonG[iZllCand[0]]->Pz());
        double angleLL = leptonG[iZllCand[1]]->Angle(a)*180/TMath::Pi();
        hDVar[iHist][73]->Fill(TMath::Min(massCloseZll,199.99),weight);
	hDVar[iHist][74]->Fill(TMath::Min(PMET,199.99),weight);
	hDVar[iHist][75]->Fill(TMath::Min((*leptonG[indexW[0]]+*leptonG[indexW[1]]).M(),199.99),weight);
        hDVar[iHist][76]->Fill(TMath::Min(leptonG[iZllCand[0]]->E()+leptonG[iZllCand[1]]->E(),199.999),weight);
	hDVar[iHist][77]->Fill(angleLL,weight); 	
	hDVar[iHist][78]->Fill(phiLL,weight);
        if(PMET > 15.0 && angleLL > 120 &&
	   (*leptonG[indexW[0]]+*leptonG[indexW[1]]).M() < 85){
	  if(abs(lType[indexW[0]]) != abs(lType[indexW[1]])){
            hDVar[iHist][79]->Fill(TMath::Min(mRecoil[0],159.999),weight);
	    if(mRecoil[0] > 120 && mRecoil[0] < 140){
	      nClass[9][iClass] += weight;
	    }
	  }
	  else {
            hDVar[iHist][80]->Fill(TMath::Min(mRecoil[0],159.999),weight);
	    if(mRecoil[0] > 120 && mRecoil[0] < 140){
	      nClass[10][iClass] += weight;
	    }
	  }	  
	  if(ifile == 0) hDSig[8]->Fill(processID,weight);
          hDVar[iHist][94]->Fill(angleLMET,weight);
        }
      }

      // ZH->nnqqln
      if(nLeptons == 1 && nJets == 2 && chargeMulti > 10 && etaJetMax < 2.5 && jetBtagMax[0] < 4.0){
	double wLQ[4] = {theMET.Px()+leptonG[indexW[0]]->Px(),theMET.Py()+leptonG[indexW[0]]->Py(),
	                 theMET.Pz()+leptonG[indexW[0]]->Pz(),theMET.Mag()+leptonG[indexW[0]]->P()};
        double massLQ = wLQ[3]*wLQ[3]-wLQ[0]*wLQ[0]-wLQ[1]*wLQ[1]-wLQ[2]*wLQ[2];
	if(massLQ > 0 && leptonG[indexW[0]]->Pt() < 40) massLQ = sqrt(massLQ);
        hDVar[iHist][81]->Fill(TMath::Min(MissingMass,199.99),weight);
        hDVar[iHist][82]->Fill(TMath::Min(chargeMulti,99.499),weight);
	hDVar[iHist][83]->Fill(TMath::Min(PMET,199.99),weight);
        hDVar[iHist][84]->Fill(jetBtagMax[0],weight);
        hDVar[iHist][85]->Fill(jetBtagMax[1],weight);
        hDVar[iHist][86]->Fill(anglesQQ[0],weight);
        hDVar[iHist][87]->Fill(angleQL,weight);
    	hDVar[iHist][88]->Fill(TMath::Min(leptonG[indexW[0]]->Pt(),199.99),weight);
	if(anglesQQ[0] > 115 && leptonG[indexW[0]]->Pt() < 40 && PMET > 10 &&
	   massLQ < 70 && massCloseWqqNoZRemoval < 100 && MissingMass > 110){
	  hDVar[iHist][89]->Fill(TMath::Min(mRecoil[2],159.99),weight);
	  if(ifile == 0) hDSig[9]->Fill(processID,weight);
          hDVar[iHist][95]->Fill(angleLMET,weight);
	  if(mRecoil[2] > 120){
	    nClass[11][iClass] += weight;
	  }
	}
      }

    } // end tree loop
    file->Close();
  } // end chain process

  for(int i=0; i<12; i++){
    Double_t nHbck = nClass[i][2] + nClass[i][3] + nClass[i][4] + nClass[i][5] + nClass[i][6] + nClass[i][7] + nClass[i][8] + nClass[i][9];
    printf("class(%10s): HWW: %7.2f HXX: %7.2f BCK: %9.2f -> precision: %5.3f\n",ClassName[i].Data(),nClass[i][1],nHbck,nClass[i][0],sqrt(nClass[i][1]+nHbck+nClass[i][0])/nClass[i][1]);
  }
  printf("-----------------------------------------------------\n");
  for(int i=0; i<12; i++){
    printf("class(%10s): Hbb: %7.2f Hcc: %7.2f Hss: %7.2f Htt: %7.2f Hmm: %7.2f HGG: %7.2f Hgg: %7.2f Hzz: %7.2f\n",ClassName[i].Data(),
           nClass[i][2],nClass[i][3],nClass[i][4],nClass[i][5],nClass[i][6],nClass[i][7],nClass[i][8],nClass[i][9]);
  }


  TFile* outFilePlotsNote = new TFile("output.root","recreate");
  outFilePlotsNote->cd();
  for(UInt_t i=0; i<10; i++){
    hDSig[i]->Write();
  }
  for(UInt_t i=0; i<infilenamev.size()+8; i++){
    for(UInt_t j=0; j<=99; j++){
      hDVar[i][j]->Write();
    }
  }
  outFilePlotsNote->Close();

}
double DeltaPhi(double phi1, double phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].
  double dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}
