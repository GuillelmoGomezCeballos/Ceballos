#!/bin/sh

export TAG=$1;
export ERA=$2;
export VERSION=$3;

echo "**************************************"
echo "TAG: $TAG ERA: $ERA VERSION: $VERSION:"
echo "**************************************"

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-mueg-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-mueg-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/$TAG/$ERA-mueg-$VERSION/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-dmu-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-dmu-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/$TAG/$ERA-dmu-$VERSION/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-smu-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-smu-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/$TAG/$ERA-smu-$VERSION/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-del-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-del-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/$TAG/$ERA-del-$VERSION/Filesets|awk '{printf("%3d\n",$1)}';

ls -l               test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-sel-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
grep ysis::Terminat test/afs/cern.ch/user/c/ceballos/scratch0/catalog/cern/filefi/$TAG/$ERA-sel-$VERSION/*out|wc|awk '{printf("%3d ",$1)}';
wc ~/scratch0/catalog/cern/filefi/$TAG/$ERA-sel-$VERSION/Filesets|awk '{printf("%3d\n",$1)}';
