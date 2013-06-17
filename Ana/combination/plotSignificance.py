#! /usr/bin/env python

from ROOT import gROOT, TLegend, TCanvas, TFile, TGraph, TMultiGraph, TLatex
import sys, os, math, argparse
import rootlogonTDR
from array import array

parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),description="analysis")
parser.add_argument('-t',dest='tag',action='store',required=True,help='tag')
parser.add_argument('-d',dest='dir',action='store',required=True,help='dir')
args=parser.parse_args()
#if args.tag not in ['8TeV_ofshape','8TeV_ofshape0','8TeV_ofshape1','8TeV_sf_cut','8TeV_bdt','8TeV_cut','8TeV_vbf_cut','8TeV_vbf_bdt','7p8TeV','wh_3l_cut','wh_3l_bdt']: args.tag='8TeV'

masses = [110,115,120,125,130,135,140,150]
masses = [110,115,120,125,130,135,140,150,160,170,180,190,200]
masses = [110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600]

significances = {}
graphs = {}
multiGraph = TMultiGraph()
legend = TLegend(0.7,0.7,0.9,0.9)
legend.SetFillColor(0)

for limitType in ['Expected', 'Observed']:
    significances[limitType] = []
    graphs[limitType] = TGraph

    for mass in masses:
        print             args.dir+'/'+str(mass)+'/higgsCombine'+limitType+args.tag+'.ProfileLikelihood.mH'+str(mass)+'.root'
        inputFile = TFile(args.dir+'/'+str(mass)+'/higgsCombine'+limitType+args.tag+'.ProfileLikelihood.mH'+str(mass)+'.root')

        tree = inputFile.Get('limit')
        tree.GetEntry(0)
        significances[limitType].append(tree.GetLeaf('limit').GetValue())
    print significances[limitType]
    graphs[limitType] = TGraph(len(masses), array('f',masses), array('f',significances[limitType]))
    graphs[limitType].SetLineWidth(2)
    if limitType=='Expected': graphs[limitType].SetLineStyle(2)
    else: graphs[limitType].SetLineStyle(1)
    multiGraph.Add(graphs[limitType])
    legend.AddEntry(graphs[limitType],limitType,'l')
    
cCanvas = TCanvas('significance')
multiGraph.Draw('aL')
multiGraph.GetXaxis().SetTitle('M_{H} (GeV)')
multiGraph.GetYaxis().SetTitle('significance')

xCorner=140
yCorner=7.5
text = [TLatex(xCorner,yCorner,'CMS Preliminary')]
if args.tag.find('7p8TeV')>-1: 
    text.append(TLatex(xCorner,yCorner-1,'#sqrt{s} = 7 TeV, L = 4.9 fb^{-1}'))
    text.append(TLatex(xCorner,yCorner-2,'#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}'))
elif args.tag.find('8')>-1: 
    text.append(TLatex(xCorner,yCorner-1,'#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}'))
elif args.tag.find('7')>-1: 
    text.append(TLatex(xCorner,yCorner-1,'#sqrt{s} = 7 TeV, L = 4.9 fb^{-1}'))
elif args.tag.find('3l')>-1: 
    text.append(TLatex(xCorner,yCorner-0.1,'#sqrt{s} = 7 TeV, L = 4.9 fb^{-1}'))
    text.append(TLatex(xCorner,yCorner-0.2,'#sqrt{s} = 8 TeV, L = 19.5 fb^{-1}'))
else: lumi='L=xxx'

legend.Draw()
for t in text:
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.Draw('same')
cCanvas.Update()

#cCanvas.SaveAs('significance'+args.tag+'.pdf')
#cCanvas.SaveAs('significance'+args.tag+'.eps')
raw_input('type whatever to quit')
