#! /usr/bin/env python

from ROOT import gROOT, TLegend, TCanvas, TFile, TGraph, TMultiGraph
import sys, os, math, argparse
import rootlogonTDR
from array import array

from drawFunction import *

parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),description="analysis")
parser.add_argument('-t',dest='tag',action='store',required=True,help='tag')
parser.add_argument('-d',dest='dir',action='store',required=True,help='dir')
args=parser.parse_args()
#if args.tag not in ['8TeV_ofshape','8TeV_sf_cut','8TeV_bdt','8TeV_cut','8TeV_vbf_cut','8TeV_vbf_bdt','7p8TeV','wh_3l_cut','wh_3l_bdt']: args.tag='8TeV'

masses = [110,115,120,125,130,135,140,150]
masses = [110,115,120,125,130,135,140,150,160,170,180,190,200]
masses = [110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600]
#masses = [105,115,125,135,145]

limits = {}

for mass in masses:
    limits[mass] = []
    print             args.dir+'/'+str(mass)+'/higgsCombinehWW'+args.tag+'.MaxLikelihoodFit.mH'+str(mass)+'.root'
    inputFile = TFile(args.dir+'/'+str(mass)+'/higgsCombinehWW'+args.tag+'.MaxLikelihoodFit.mH'+str(mass)+'.root')

    tree = inputFile.Get('limit')
    # expected
    tree.GetEntry(0)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    #-s
    tree.GetEntry(1)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    #+s
    tree.GetEntry(2)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    # observed
    tree.GetEntry(3)
    limits[mass].append(tree.GetLeaf('limit').GetValue())

#print limits

drawBestFit(limits, args.tag)

