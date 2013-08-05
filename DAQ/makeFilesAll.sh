#!/bin/sh

rm -rf unmerged/* merged/* done/*

./makeFiles.py --streamName=STREAMA --sizePerFile=1 &

./makeFiles.py --streamName=STREAMB --sizePerFile=5 &

./makeFiles.py --streamName=STREAMC --sizePerFile=20 &
