#!/usr/bin/env python

import os, time, sys, getopt
import random

"""
Do actual files
"""
def doFiles(RUNNumber):
   random.seed(int(seeds))
   theNLoop = 1
   nInput = 0
   nOutput = 0
   LSNumber = 0

   start = time.time()
   while (float(timeEnd) < 0.0 or float(time.time() - start) < float(timeEnd)):
     time.sleep (float(rate))

     # just in case we need more than one seed
     numberOfSeedsNeeded = 1
     seedsRND = []
     for i in range(0, numberOfSeedsNeeded):
       seedsRND.append(random.randint(0,999999))

     nInput += 100
     nOutput += 9
     fileName =  "%s/Data.%d.LS%d.%s.%d.std.dat" % (path_to_make,RUNNumber,LSNumber,streamName,seedsRND[0])
     thefile = open(fileName, 'w')
     thefile.write('0' * 1024 * 1024 * int(sizePerFile))
     thefile.write("\n")
     #msg  = "A B C D E F G H I - 0 - %d\n" % theNLoop
     #for j in range(0, 50000000):
     #   msg += "%d " % j
     #   if(j%1000 == 0):
     #	   msg += "\n"
     #thefile.write(msg)
     thefile.close()

     fileName =  "%s/Data.%d.LS%d.%s.%d.mon.dat" % (path_to_make,RUNNumber,LSNumber,streamName,seedsRND[0])
     thefile = open(fileName, 'w')
     thefile.write('0' * 1024 * 1024 * 1)
     thefile.write("\n")
     #msg = "A B C D E F G H I - 1 - %d\n" % theNLoop
     #for j in range(0, 2000000):
     #   msg += "%d " % j
     #   if(j%1000 == 0):
     #	   msg += "\n"
     #thefile.write(msg)
     thefile.close()

     fileName =  "%s/Data.%d.LS%d.%s.%d.met.dat" % (path_to_make,RUNNumber,LSNumber,streamName,seedsRND[0])
     thefile = open(fileName, 'w')
     thefile.write("100 9\n")
     thefile.close()

     if(theNLoop%30 == 0):
	fileName =  "%s/Data.%d.LS%d.%s.0000.eol.dat" % (path_to_make,RUNNumber,LSNumber,streamName)
	thefile = open(fileName, 'w')
	thefile.write("%f %f\n" % (nInput,nOutput))
	thefile.close()
	nInput = 0
	nOutput = 0
	LSNumber += 1

     if(theNLoop%1000 == 0):
	nInput = 0
	nOutput = 0
	LSNumber += 1
	RUNNumber += 1

     theNLoop += 1

"""
Main
"""
valid = ['path_to_make=', 'rate=', 'seeds=', 'timeEnd=', 
         'run=', 'streamName=', 'sizePerFile=', 
         'help']

usage =  "Usage: listdir.py --path_to_make=<path to write files>\n"
usage += "                  --rate=<rate_to_create_files (in seconds)>\n"
usage += "                  --seeds=<initial seed>\n"
usage += "                  --timeEnd=<time until stop (in seconds, 0 will never stop)>\n"
usage += "                  --streamName=<stream name, (STREAMA by default)>\n"
usage += "                  --sizePerFile=<in MB (50 by default, must be an integer)>\n"
usage += "                  --run=<initial run>\n"

try:
   opts, args = getopt.getopt(sys.argv[1:], "", valid)
except getopt.GetoptError, ex:
   print usage
   print str(ex)
   sys.exit(1)

path_to_make = "unmerged"
rate         = 0.0
seeds        = 999
timeEnd      = -1
RUNNumber    = 100
streamName   = "STREAMA"
sizePerFile  = 50

for opt, arg in opts:
   if opt == "--help":
      print usage
      sys.exit(1)
   if opt == "--path_to_make":
      path_to_make = arg
   if opt == "--rate":
      rate = arg
   if opt == "--seeds":
      seeds = arg
   if opt == "--timeEnd":
      timeEnd = arg
   if opt == "--run":
      RUNNumber = arg
   if opt == "--streamName":
      streamName = arg
   if opt == "--sizePerFile":
      sizePerFile = arg

if not os.path.exists(path_to_make):
   os.mkdir(path_to_make)

doFiles(int(RUNNumber))
