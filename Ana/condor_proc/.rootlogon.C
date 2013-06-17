{
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }

  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  loadLibraries("libMitPhysics*.so");
  loadLibraries("libMitHiggs*.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libAnaSelMods.so");
}
