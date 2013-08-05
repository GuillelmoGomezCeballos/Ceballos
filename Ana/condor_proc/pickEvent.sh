stager_qry -M /castor/cern.ch/cms/store/data/Run2011B/MuEG/AOD/PromptReco-v1/000/178/803/BCF4FF19-35FA-E011-A0E6-001D09F2462D.root
stager_qry -M /castor/cern.ch/cms/store/data/Run2011B/MuEG/AOD/PromptReco-v1/000/178/920/8A502F54-E6FB-E011-88A1-BCAEC5329728.root
edmCopyPickMerge outputFile=pickevents_run_178803_VBF.root eventsToProcess=178803:449292809 inputFiles=/store/data/Run2011B/MuEG/AOD/PromptReco-v1/000/178/803/BCF4FF19-35FA-E011-A0E6-001D09F2462D.root
edmCopyPickMerge outputFile=pickevents_run_178920_VBF.root eventsToProcess=178920:118380801 inputFiles=/store/data/Run2011B/MuEG/AOD/PromptReco-v1/000/178/920/8A502F54-E6FB-E011-88A1-BCAEC5329728.root

