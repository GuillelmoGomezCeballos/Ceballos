//#include <iostream>
//#include <fstream> 
//#include <sstream>
//#include "TFile.h"
//#include "TString.h"
//#include "RooArgSet.h"
//#include "RooRealVar.h"

using namespace std;

void printtableline(TString process, float Npre, float Npost) {
    cout.width(15); cout << process << " & ";
    cout.width(10); cout << Form("%5.1f", Npre) << " & ";
    cout.width(10); cout << Form("%5.1f", Npost) << " & ";
    cout.width(10); cout << Form("%5.1f", (Npost - Npre)) << " & ";
    cout.width(10); if(Npre > 0) cout << Form("%5.1f", (Npost - Npre)/Npre*100);
                    else         cout << Form("%5.1f", 0.0);
    cout.width(10); cout << "\\\\" << endl;
}
//main
void printNorm_combine_ZH(TString mlfit="mlfit.root", TString card="cards/test/125/hwwof_0j_shape_8TeV.txt",
                          bool printsignal = true) {

    // ---------------------
    //      PRE-FIT 
    // ---------------------

    // get prefit normalization from a card   
    string line;
    ifstream incard(card.Data()); 
    
    char* tmp   =   "";
    float pre_ZH     =   0.;	 float post_ZH     =	0.;
    float pre_Zjets  =   0.;	 float post_Zjets  =	0.;
    float pre_VVV    =   0.;	 float post_VVV    =	0.;
    float pre_ZZ     =   0.;	 float post_ZZ     =	0.;
    float pre_WZ     =   0.;	 float post_WZ     =	0.;
    float pre_EM     =   0.;	 float post_EM     =	0.;
    float pre_Wjets  =   0.;	 float post_Wjets  =	0.;

    if (incard.is_open())
    {
        while ( incard.good() )
        { 
            // get a line from input file
            getline (incard,line);

            // skip the line if it does not contain "rate" 
            if( line.find("rate")==string::npos ) continue;
              
            stringstream stream(line);
            stream >> tmp;
            stream >> pre_ZH;
            stream >> pre_Zjets;
            stream >> pre_VVV;
            stream >> pre_WZ;
            stream >> pre_ZZ;
            stream >> pre_EM; 
            stream >> pre_Wjets;

            if( !incard.good() ) continue;
        }
    } else { cout << "problem reading, exiting..." << endl; return;}
    incard.close(); 

    cout << "have read prefit yields..." << endl;

    // ---------------------
    //      POST-FIT 
    // ---------------------
    TFile *postFile = TFile::Open(mlfit.Data(), "READ"); 
    postFile->ls(); 
    RooArgSet *norm_fit_s = (RooArgSet*)postFile->Get("norm_fit_s");
    norm_fit_s->Print("v"); // for debugging  
    RooRealVar *r; 
    
    if(pre_ZH)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_ZH");     post_ZH    = (float)r->getVal(); } 
    if(pre_Zjets)   { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_Zjets");        post_Zjets = (float)r->getVal(); } 
    if(pre_VVV)     { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_VVV");    post_VVV   = (float)r->getVal(); } 
    if(pre_WZ)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_WZ");     post_WZ    = (float)r->getVal(); } 
    if(pre_ZZ)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_ZZ");     post_ZZ    = (float)r->getVal(); } 
    if(pre_EM)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_EM");     post_EM    = (float)r->getVal(); } 
    if(pre_Wjets)   { r =  (RooRealVar*) norm_fit_s->find("n_exp_binch2_proc_Wjets");  post_Wjets = (float)r->getVal(); } 

    // 
    // print table 
    // 
    // -------------------------------------------------------------------------------
    //      print table 
    // -------------------------------------------------------------------------------
    cout << "\\begin{table}[ht!]" << endl;
    cout << "\\begin{center}" << endl;
    cout << "\\begin{tabular}{c|cc|cc}" << endl;
    cout << "\\hline" << endl;   cout << "\\hline" << endl;
    cout << "        Process &    N(prefit) &   N(postfit) & Difference(raw) &  Difference(\\%)  \\\\  " << endl;
    cout << "\\hline" << endl;   cout << "\\hline" << endl;
    if(printsignal) {
        printtableline("ZH",     pre_ZH,     post_ZH);
        cout << "\\hline" << endl;
    }
    printtableline("Zjets",   pre_Zjets,   post_Zjets);
    printtableline("VVV",     pre_VVV,     post_VVV);
    printtableline("WZ",      pre_WZ,      post_WZ);
    printtableline("ZZ",      pre_ZZ,      post_ZZ);
    printtableline("WW/top",  pre_EM,      post_EM);
    printtableline("W+jets",  pre_Wjets,   post_Wjets);
    cout << "\\hline" << endl;
    cout << "\\caption{ }" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\end{table}" << endl;
}
