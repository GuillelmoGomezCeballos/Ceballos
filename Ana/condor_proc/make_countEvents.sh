#!/bin/sh

cp countEvents_begin.txt countEvents.C;
awk -f $HOME/condor/make_countEvent.awk submitter_all.csh >> countEvents.C;
cat countEvents_end.txt >> countEvents.C;
export NFILES=`wc -l submitter_all.csh|awk '{print$1}'`;
sed -ie "s/99999/${NFILES}/" countEvents.C;
rm -f countEvents.Ce;
