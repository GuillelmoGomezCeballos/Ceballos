#!/bin/tcsh
echo "Arguments: $*"

echo "*********** working node ***********";
/bin/hostname;
echo "*********** id ***********";
/usr/bin/id
echo "*********** env ***********";
/usr/bin/env

echo " "; echo "Initialize CMSSW"; echo " "
pwd
set pwd=`pwd`
#source $HOME/EVAL64 4_2_8_patch4
source $HOME/EVAL65 5_2_3_patch3
cd     $HOME/condor
setenv KRB5CCNAME $HOME/.krb5/ticket

set logFile=/tmp/mylog_$$.log

echo " "; echo "Starting root now"; echo " "
echo "root -l -b -q .rootlogon.C $*";
#root -l -b -q .rootlogon.C $*;
set STAGE_SVCCLASS=cmsprod;
root -b -l -q .rootlogon.C $1\(\"$2\",\"$3\",\"$4\",\"$5\",$6,$7\);

#set status=`echo $?`
#set error=`grep Error $logFile | wc -l`
#set zip=`grep R__unzip: $logFile | wc -l`
  
#echo "Status: $status  Errors: $error  R__Unzip: $zip"
#if ( $status == 0 && $error == 0 && $zip == 0 ) then
#  echo "DECISION - looks healthy."
#  echo " "
#else
#  echo "DECISION - looks corrupted."
#endif

exit 0;
