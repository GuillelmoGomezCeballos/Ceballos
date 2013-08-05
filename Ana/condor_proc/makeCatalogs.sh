#!/bin/sh

echo '***** catalogs: 025 *****';
grep 025 submitter_all.csh > temp_submitter_all.csh;
awk '{print"rm -rf ~/scratch0/catalog/"$4"/"$5";mkdir -p ~/scratch0/catalog/"$4";scp -r dtmit05:/home/mitprod/catalog/"$4"/"$5" ~/scratch0/catalog/"$4"/"$5}' temp_submitter_all.csh > temp;
source temp;
rm -f temp;

awk '{print"wc ~/scratch0/catalog/"$4"/"$5"/Files;rfdir /castor/cern.ch/user/p/paus/filefi/025/"$5" | grep root | wc;"}' temp_submitter_all.csh >> temp;
source temp | awk '{printf("%d ",$1); if(NR%2==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD "$1" "$2};'
rm -f temp;
rm -f temp_submitter_all.csh;

echo '***** catalogs: 028 *****';
grep 028 submitter_all.csh > temp_submitter_all.csh;
awk '{print"rm -rf ~/scratch0/catalog/"$4"/"$5";mkdir -p ~/scratch0/catalog/"$4";scp -r dtmit05:/home/mitprod/catalog/"$4"/"$5" ~/scratch0/catalog/"$4"/"$5}' temp_submitter_all.csh > temp;
source temp;
rm -f temp;

awk '{print"wc ~/scratch0/catalog/"$4"/"$5"/Files;rfdir /castor/cern.ch/user/p/paus/filefi/028/"$5" | grep root | wc;"}' temp_submitter_all.csh >> temp;
source temp | awk '{printf("%d ",$1); if(NR%2==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD "$1" "$2};'
rm -f temp;
rm -f temp_submitter_all.csh;

echo '***** catalogs: 029 *****';
grep 029 submitter_all.csh > temp_submitter_all.csh;
awk '{print"rm -rf ~/scratch0/catalog/"$4"/"$5";mkdir -p ~/scratch0/catalog/"$4";scp -r dtmit05:/home/mitprod/catalog/"$4"/"$5" ~/scratch0/catalog/"$4"/"$5}' temp_submitter_all.csh > temp;
source temp;
rm -f temp;

awk '{print"wc ~/scratch0/catalog/"$4"/"$5"/Files;rfdir /castor/cern.ch/user/p/paus/filefi/029/"$5" | grep root | wc;"}' temp_submitter_all.csh >> temp;
source temp | awk '{printf("%d ",$1); if(NR%2==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD "$1" "$2};'
rm -f temp;
rm -f temp_submitter_all.csh;

echo '***** catalogs: 030 *****';
grep 030 submitter_all.csh > temp_submitter_all.csh;
awk '{print"rm -rf ~/scratch0/catalog/"$4"/"$5";mkdir -p ~/scratch0/catalog/"$4";scp -r dtmit05:/home/mitprod/catalog/"$4"/"$5" ~/scratch0/catalog/"$4"/"$5}' temp_submitter_all.csh > temp;
source temp;
rm -f temp;

awk '{print"wc ~/scratch0/catalog/"$4"/"$5"/Files;rfdir /castor/cern.ch/user/p/paus/filefi/030/"$5" | grep root | wc;"}' temp_submitter_all.csh >> temp;
source temp | awk '{printf("%d ",$1); if(NR%2==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD "$1" "$2};'
rm -f temp;
rm -f temp_submitter_all.csh;

echo '***** catalogs: 031 *****';
grep 031 submitter_all.csh > temp_submitter_all.csh;
awk '{print"rm -rf ~/scratch0/catalog/"$4"/"$5";mkdir -p ~/scratch0/catalog/"$4";scp -r dtmit05:/home/mitprod/catalog/"$4"/"$5" ~/scratch0/catalog/"$4"/"$5}' temp_submitter_all.csh > temp;
source temp;
rm -f temp;

awk '{print"wc ~/scratch0/catalog/"$4"/"$5"/Files;rfdir /castor/cern.ch/user/p/paus/filefi/031/"$5" | grep root | wc;"}' temp_submitter_all.csh >> temp;
source temp | awk '{printf("%d ",$1); if(NR%2==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD "$1" "$2};'
rm -f temp;
rm -f temp_submitter_all.csh;

rm -rf ~/scratch0/catalog/*/*/*/*/old;
