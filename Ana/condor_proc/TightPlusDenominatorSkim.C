// root -l $CMSSW_BASE/src/MitPhysics/macros/skimming/TightPlusDenominatorSkim.C+\(\"0000\",\"r10b-el-pr-v2\",\"cern/filefi/014/\",\"/home/mitprod/catalog/\",10000\)
// $Id: TightPlusDenominatorSkim.C,v 1.1 2012/10/05 09:06:41 ceballos Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod/interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void TightPlusDenominatorSkim(const char *catalogDir   = "/home/ceballos/catalog",
                  const char *book	   = "/cern/filler/014a",
                  const char *dataset	   = "p10-wjets-mg-v26",
                  const char *fileset	   = "0003",
                  int nsel = 0, int nEvents = 1000)
{
  TString skimName("singlelepton");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muId = new MuonIDMod;  
  muId->SetPtMin     (17.0);
  muId->SetApplyD0Cut(kFALSE);
  muId->SetClassType ("Global");
  muId->SetIDType    ("Minimal");
  muId->SetIsoType   ("TrackCaloSliding");

  ElectronIDMod *electronNoId = new ElectronIDMod;
  electronNoId->SetInputName                 (Names::gkElectronBrn);
  electronNoId->SetPtMin                     (10.0);
  electronNoId->SetApplyConversionFilterType1(kFALSE);
  electronNoId->SetApplyConversionFilterType2(kFALSE);
  electronNoId->SetChargeFilter              (kFALSE);
  electronNoId->SetApplySpikeRemoval         (kFALSE);
  electronNoId->SetApplyD0Cut                (kFALSE);
  electronNoId->SetNExpectedHitsInnerCut     (999);
  electronNoId->SetIDType                    ("NoId");
  electronNoId->SetIsoType                   ("NoIso");

  ElectronIDMod       *electronIDWP80IdIso       = new ElectronIDMod;
  electronIDWP80IdIso->SetPtMin(0.0);
  electronIDWP80IdIso->SetEtaMax(999.);
  electronIDWP80IdIso->SetEtMin(17.0);
  electronIDWP80IdIso->SetIDType(TString("VBTFWorkingPoint80Id"));
  electronIDWP80IdIso->SetIsoType(TString("VBTFWorkingPoint80Iso"));
  electronIDWP80IdIso->SetNExpectedHitsInnerCut(0);
  electronIDWP80IdIso->SetApplyConversionFilterType1(kFALSE);
  electronIDWP80IdIso->SetApplyConversionFilterType2(kTRUE);
  electronIDWP80IdIso->SetApplyD0Cut(kFALSE);
  electronIDWP80IdIso->SetChargeFilter(kFALSE);
  electronIDWP80IdIso->SetApplyTriggerMatching(kFALSE);
  electronIDWP80IdIso->SetApplyEcalFiducial(kTRUE);
  electronIDWP80IdIso->SetApplyEcalSeeded(kTRUE);
  electronIDWP80IdIso->SetApplyCombinedIso(kTRUE);
  electronIDWP80IdIso->SetApplySpikeRemoval(kTRUE);
  electronIDWP80IdIso->SetOutputName("GoodElectronsWP80IdIso");


  MergeLeptonsMod *mergedTight = new MergeLeptonsMod;
  mergedTight->SetMuonsName    (muId->GetOutputName());
  mergedTight->SetElectronsName(electronIDWP80IdIso->GetOutputName());
  mergedTight->SetOutputName("mergedTightLeptons");

  MergeLeptonsMod *mergedLoose = new MergeLeptonsMod;
  mergedLoose->SetMuonsName    (muId->GetOutputName());
  mergedLoose->SetElectronsName(electronNoId->GetOutputName());
  mergedLoose->SetOutputName("mergedLooseLeptons");

  GenericSelMod<mithep::Particle> *selModTight = new GenericSelMod<mithep::Particle>;
  selModTight->SetPtMin(0.0);
  selModTight->SetMinCounts(1);
  selModTight->SetColName(mergedTight->GetOutputName());

  GenericSelMod<mithep::Particle> *selModLoose = new GenericSelMod<mithep::Particle>;
  selModLoose->SetPtMin(0.0);
  selModLoose->SetMinCounts(2);
  selModLoose->SetColName(mergedLoose->GetOutputName());



  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");

  TString rootFile = "";
  rootFile += TString("skims/") + dataset + TString("/");
  rootFile += skimName;
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  outMod->SetFileName(rootFile);
  outMod->SetPathName(".");

  //------------------------------------------------------------------------------------------------
  // set up analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->AddSuperModule(muId);
  muId->Add(electronIDWP80IdIso);
  electronIDWP80IdIso->Add(electronNoId);
  electronNoId->Add(mergedTight);
  mergedTight->Add(mergedLoose);
  mergedLoose->Add(selModTight);
  selModTight->Add(selModLoose);
  selModLoose->Add(outMod); 
  if (nEvents>0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  Catalog *c = new Catalog(catalogDir);
  Dataset *d = c->FindDataset(book,dataset,fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(kFALSE);
}
