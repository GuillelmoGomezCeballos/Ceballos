#!/bin/sh

./checkFiles.sh | awk '(NR%3!=1){printf("%d ",$1); if(NR%3==0)printf("\n");}' | awk '{if ($1-$2==0) print "GOOD";else print NR" BADDDDDDDDDDD: ",$1,$2};'
