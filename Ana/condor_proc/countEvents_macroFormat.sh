#!/bin/bash

root -l -q -b countEvents.C >& log;
grep ISR log | awk '{split($1,a,"old/");split(a[2],b,".root:");print"#define EVENTS_"b[1]"   "$2}';
rm log;
