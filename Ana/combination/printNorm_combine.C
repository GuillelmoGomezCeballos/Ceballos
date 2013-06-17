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
    cout.width(10); cout << Form("%5.1f", (Npre ? Npost : 0.0)) << " & ";
    cout.width(10); cout << Form("%5.1f", (Npre ? (Npost - Npre) : 0.0)) << " & ";
    //cout.width(10); cout << Form("%5.1f", (Npre ? (Npost - Npre)/Npre*100 : 0.0));
    cout.width(10); cout << Form("%5.3f", (Npre ? Npost/Npre : 0.0));
    cout.width(10); cout << "\\\\" << endl;
}
//main
void printNorm_combine(TString mlfit="mlfit.root", TString card="cards/test/125/hwwof_0j_shape_8TeV.txt") {

    // ---------------------
    //      PRE-FIT 
    // ---------------------

    // get prefit normalization from a card   
    string line;
    ifstream incard(card.Data()); 
    
    char* tmp   =   "";
    float pre_ZH      =   0.;	 float post_ZH      =	0.;
    float pre_WH      =   0.;	 float post_WH      =	0.;
    float pre_qqH     =   0.;	 float post_qqH     =	0.;
    float pre_ggH     =   0.;	 float post_ggH     =	0.;
    float pre_qqWW    =   0.;	 float post_qqWW    =	0.;
    float pre_ggWW    =   0.;	 float post_ggWW    =	0.;
    float pre_VV      =   0.;	 float post_VV      =	0.;
    float pre_Top     =   0.;	 float post_Top     =	0.;
    float pre_Zjets   =   0.;	 float post_Zjets   =	0.;
    float pre_WjetsE  =   0.;	 float post_WjetsE  =	0.;
    float pre_WjetsM  =   0.;	 float post_WjetsM  =	0.;
    float pre_Wgamma  =   0.;	 float post_Wgamma  =	0.;
    float pre_Wg3l    =   0.;	 float post_Wg3l    =	0.;
    float pre_Ztt     =   0.;	 float post_Ztt     =	0.;
    float pre_ZH_SM   =   0.;	 float post_ZH_SM   =	0.;
    float pre_WH_SM   =   0.;	 float post_WH_SM   =	0.;
    float pre_qqH_SM  =   0.;	 float post_qqH_SM  =	0.;
    float pre_ggH_SM  =   0.;	 float post_ggH_SM  =	0.;
    float pre_ggH_ALT =   0.;	 float post_ggH_ALT =	0.;

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
            stream >> pre_ggH_ALT;
            stream >> pre_ZH;
            stream >> pre_WH;
            stream >> pre_qqH;
            stream >> pre_ggH;
            stream >> pre_qqWW;
            stream >> pre_ggWW;
            stream >> pre_VV;
            stream >> pre_Top;
            stream >> pre_Zjets;
            stream >> pre_WjetsE;
            stream >> pre_Wgamma;
            stream >> pre_Wg3l;
            stream >> pre_Ztt;
            stream >> pre_WjetsM; 
            stream >> pre_ZH_SM;
            stream >> pre_WH_SM;
            stream >> pre_qqH_SM;
            stream >> pre_ggH_SM;

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
    
    //if(pre_ggH_ALT)  { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_ggH_ALT"); post_ggH_ALT =	(float)r->getVal(); } 
    if(pre_ZH)       { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_ZH");      post_ZH      =	(float)r->getVal(); } 
    if(pre_WH)       { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_WH");      post_WH      =	(float)r->getVal(); } 
    if(pre_qqH)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_qqH");     post_qqH     =	(float)r->getVal(); } 
    if(pre_ggH)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_ggH");     post_ggH     =	(float)r->getVal(); } 
    if(pre_qqWW)     { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_qqWW");    post_qqWW    =	(float)r->getVal(); } 
    if(pre_ggWW)     { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_ggWW");    post_ggWW    =	(float)r->getVal(); } 
    if(pre_VV)       { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_VV");      post_VV      =	(float)r->getVal(); } 
    if(pre_Top)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_Top");     post_Top     =	(float)r->getVal(); } 
    if(pre_Zjets)    { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_Zjets");   post_Zjets   =	(float)r->getVal(); } 
    if(pre_WjetsE)   { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_WjetsE");  post_WjetsE  =	(float)r->getVal(); } 
    if(pre_WjetsM)   { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_WjetsM");  post_WjetsM  =	(float)r->getVal(); } 
    if(pre_Wgamma)   { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_Wgamma");  post_Wgamma  =	(float)r->getVal(); } 
    if(pre_Wg3l)     { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_Wg3l");    post_Wg3l    =	(float)r->getVal(); } 
    if(pre_Ztt)      { r =  (RooRealVar*) norm_fit_s->find("n_exp_final_binch1_proc_Ztt");     post_Ztt     =	(float)r->getVal(); } 
    
    
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
    //if(printsignal) {
        printtableline("ZH",     pre_ZH,     post_ZH);
        printtableline("WH",     pre_WH,     post_WH);
        printtableline("qqH",    pre_qqH,    post_qqH);
	if(post_ggH > 0)
        printtableline("ggH",    pre_ggH,    post_ggH);
	if(post_ggH_ALT > 0)
        printtableline("ggH_ALT",pre_ggH_ALT,post_ggH_ALT);
        cout << "\\hline" << endl;
   // }
    printtableline("qqWW",   pre_qqWW,   post_qqWW);
    printtableline("ggWW",   pre_ggWW,   post_ggWW);
    cout << "\\hline" << endl;
    printtableline("VV",     pre_VV,     post_VV);
    cout << "\\hline" << endl;
    printtableline("Top",    pre_Top,    post_Top);
    //cout << "\\hline" << endl;
    printtableline("Zjets",  pre_Zjets,  post_Zjets);
    cout << "\\hline" << endl;
    printtableline("Wjets($e$)", pre_WjetsE, post_WjetsE);
    printtableline("Wjets($\\mu$)", pre_WjetsM, post_WjetsM);
    cout << "\\hline" << endl;
    printtableline("W$\\gamma$", pre_Wgamma, post_Wgamma);
    printtableline("W$\\gamma$*",   pre_Wg3l,   post_Wg3l);
    cout << "\\hline" << endl;
    printtableline("Ztt",    pre_Ztt,    post_Ztt);
    cout << "\\hline" << endl;   cout << "\\hline" << endl;
    cout << "\\end{tabular}" << endl;
    //cout << "\\caption{ }" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\end{table}" << endl;
}
