#!/bin/sh

awk '{split($4,a,"cern");print"~ceballos/public/bin/staging.csh /castor/cern.ch/user/p/paus"a[2]"/"$5";sleep 200;"}' submitter_all.csh > staging_all.sh;
chmod a+x staging_all.sh;
