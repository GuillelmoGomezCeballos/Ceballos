#!/bin/sh

echo "*** INIT ***";
bjobs|grep cebal > kkk1;
#awk '{print"bpeek "$1"|grep AddFil"}' kkk1 > kkk2;
awk '{print"bpeek "$1"|grep Adding"}' kkk1 > kkk2;

echo "*** step2 ***";
#source kkk2|awk '{print"stager_qry -M "$5}' > kkk3;
source kkk2|awk '{split($5,a,"root://castorcms/");print"stager_qry -M "a[2]}' > kkk3;

echo "*** step3 ***";
source kkk3|grep STAGEIN|awk '{split($1,a,"/");print"grep -v "a[10]" ~/scratch0/catalog/cern/filefi/"a[8]"/"a[9]"/Files > Files_tmp;wc Files_tmp;wc ~/scratch0/catalog/cern/filefi/"a[8]"/"a[9]"/Files;mv Files_tmp ~/scratch0/catalog/cern/filefi/"a[8]"/"a[9]"/Files;"}' > kkk4;

echo "*** step4 ***";
source kkk4;

echo "*** bjobs -all / bkill ***";
awk '{print"bjobs -all "$1"|grep Job|grep Name"}' kkk1 > kkk5;
awk '{print"bkill "$1}' kkk1 > kkk6;

echo "*** step5 ***";
source kkk5|awk '{split($5,a,"<");split(a[2],b,">,");split(b[1],c,"_");print"./submitter2.csh 0 /afs/cern.ch/user/c/ceballos/scratch0/catalog cern/filefi/029 "c[1]" "c[2]" "c[3]" 10000000 no"}' > kkk7;

sed -ie 's/stager_qry/stager_get/' kkk3;
rm kkk3e;

echo "to stage: source kkk3";
echo "to kill: source kkk6";
echo "to resubmit: source kkk7";
