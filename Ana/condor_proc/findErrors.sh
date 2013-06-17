mkdir -p $HOME/condor/test/afs
grep FindServerAndPath          $HOME/condor/test*/afs/* -r |awk '{print$1}';
grep AbortAnalysis              $HOME/condor/test*/afs/* -r |awk '{print$1}';
grep G__exception               $HOME/condor/test*/afs/* -r |awk '{print$1}';
grep segmentation               $HOME/condor/test*/afs/* -r |awk '{print$0}';
grep aborting                   $HOME/condor/test*/afs/* -r |awk '{print$0}';
grep exceeded                   $HOME/condor/test*/afs/* -r |awk '{print$0}';
grep "TXNetFile::CreateXClient" $HOME/condor/test*/afs/* -r |awk '{print$0}';
grep "SysError"                 $HOME/condor/test*/afs/* -r |awk '{print$0}';
grep "keyboard interrupt"       $HOME/condor/test*/afs/* -r |awk '{print$0}';
