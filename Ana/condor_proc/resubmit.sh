#!/bin/sh

bjobs|wc;bjobs -r|wc;
./findErrors.sh > kkk;
awk '{split($1,a,"out");split(a[1],b,"/");split(b[20],c,"_");print"./submitter2.csh 0 /afs/cern.ch/user/c/ceballos/scratch0/catalog "b[16]"/"b[17]"/"b[18]" "b[19]" "c[1]" "c[2]" "c[3]" no"}' kkk > sss;
sort -u sss > aaa;
mv aaa sss;
rm kkk;
chmod a+x sss;
wc sss;./sss;rm sss;
bjobs|wc;bjobs -r|wc;
#edit kkk;
#awk '{print$1}' kkk| sort -u > aaa;
#mv aaa kkk;
#awk '{print"rm "$1}' kkk > sss;
#edit kkk;
#awk '{print"./submitter2.csh 0 "$1" "$2" "$3" "$4" "$5" "$6" no"}' kkk > sss;
