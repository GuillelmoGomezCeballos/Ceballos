#!/bin/sh

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/025/r11a-pho-j16-v1-bp/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/025/r11a-pho-j16-v1-bp/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/025/r11a-pho-j16-v1-bp/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/025/r11b-pho-j16-v1-bp/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/025/r11b-pho-j16-v1-bp/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/025/r11b-pho-j16-v1-bp/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12a-pho-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12a-pho-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/029/r12a-pho-j22-v1/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12b-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12b-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/029/r12b-sph-j22-v1/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12c-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12c-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/029/r12c-sph-j22-v1/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12d-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/029/r12d-sph-j22-v1/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/029/r12d-sph-j22-v1/Filesets|awk '{printf("%3d\n",$1)}';
