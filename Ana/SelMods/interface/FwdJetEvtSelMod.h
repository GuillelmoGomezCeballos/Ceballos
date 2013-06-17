//-----------------------------------------------------------------------------
// $Id: FwdJetEvtSelMod.h,v 1.3 2010/10/19 22:09:42 ceballos Exp $
//
// FwdJetEvtSelMod
//
// A Module for Selecting forward tagging jets
// and produces some distributions
//
//
// Authors: ceballos
//-----------------------------------------------------------------------------

#ifndef HWWMODS_FWDJETEVTSELMOD_H
#define HWWMODS_FWDJETEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class FwdJetEvtSelMod : public BaseMod
  {
    public:
    FwdJetEvtSelMod(const char *name="FwdJetEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~FwdJetEvtSelMod() {}
      void         SetPrintDebug(bool b)            { fPrintDebug = b;	      }
      void         SetFillHist(bool b)              { fFillHist = b;	      }
      void         SetCleanJetsName (TString s)     { fCleanJetsName = s;     }
      void         SetCleanFwdJetsName (TString s)  { fCleanFwdJetsName = s;  }
      void         SetCleanNoFwdJetsName (TString s){ fCleanNoFwdJetsName = s;}
      void         SetMetName(TString s)            { fMetName = s;           }   
      void         SetMCqqHsName (TString s)        { fMCqqHsName = s;        }
      void         SetUseANN (bool b)               { fUseANN = b;            }
      void         SetUseHighestPtJets (bool b)     { fUseHighestPtJets= b;   }
      void         SetJetPtMax (double x)           { fJetPtMax = x;	      }
      void         SetJetPtMin (double x)           { fJetPtMin = x;	      }
      void         SetDeltaEtaMin (double x)        { fDeltaEtaMin = x;       }
      void         SetDiJetMassMin (double x)       { fDiJetMassMin = x;      }
      void         SetANNOutputMin (double x)       { fANNOutputMin = x;      }
      void         SetPtJetCut(double x)            { fPtJetCut = x;          }

    protected:
      bool         fPrintDebug;
      Double_t     fPtJetCut;
      bool         fFillHist;
      TString      fMetName;
      TString      fCleanJetsName; 
      TString      fCleanFwdJetsName; 
      TString      fCleanNoFwdJetsName;
      TString      fMCqqHsName;
  
      TH1D         *hDFwdJetSel[300];
      TH1D         *hDFwdJetLepSel[100];

      bool         fUseANN;
      bool         fUseHighestPtJets;
      double       fJetPtMax;
      double       fJetPtMin;
      double       fDeltaEtaMin;
      double       fDiJetMassMin;
      double       fANNOutputMin;
      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      
      double       Testmlp_qqH(double var[6]);
      double       Sigmoid(double x);

      ClassDef(FwdJetEvtSelMod,1) // TAM example analysis module
  };
}
#endif
