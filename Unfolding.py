from ROOT import gROOT
gROOT.ProcessLine( 'gSystem->Load("/Users/tuomas/OneDrive/work/032.JTAnalysis/Unfolding/Root6/RooUnfold/libRooUnfold");')
from ROOT import gRandom, TH1, TH1D, cout, TF1, TH2D, TH2
from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import TMath
from ROOT import TVector3
from ROOT import TRandom3
from ROOT import TNtuple
from ROOT import TFile
import matplotlib
import rootpy.ROOT as ROOT
# from ROOT import RooUnfoldSvd
# from ROOT import RooUnfoldTUnfold
# from ROOT import RooUnfoldIds
import rootpy
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
from rootpy.plotting import Hist,Hist2D,Hist3D
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import math
from rootpy.io import root_open
import sys
from drawing import *
from defs import *
import time
import JtUnfolder
from ctypes import c_int

def main():
  print 'Number of arguments: ', len(sys.argv), 'arguments.'
  print 'Argument list:',str(sys.argv)
  filename = sys.argv[1]
  print "Input file: {}".format(filename)
  
  jetBinBorders = [5,10,20,30,40,60,80,100,150,500]
  jetPtBins = [(a,b) for a,b in zip(jetBinBorders,jetBinBorders[1:])]
  JetPtCenter = [7.5,15,25,35,50,70,90,125,325]
  JetPtError = [2.5,5,5,5,10,10,10,25,175]
  Njets = len(jetBinBorders)-1
  Njets = 8
  f = root_open(filename, 'read')
  nR = 2
  iFinder = 1
  iMCFinder = iFinder + nR

  numberJetsMeas = [f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)).GetEntries() for iJet in range(Njets)]
  numberJetsTrue = [f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iMCFinder,iJet)).GetEntries() for iJet in range(Njets)]
  hTrackJtMeas = [f.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)] 
  hTrackJtTrue = [f.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iMCFinder,iJet)) for iJet in range(Njets)] 

  
  print("Measured Jets:")
  print(numberJetsMeas)
  print("True Jets:")
  print(numberJetsTrue)
  hTrackJtCorr2D = [f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackJtCorr2D/TrackJtCorr2DNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)]
  hTrackJtMeas2D = f.Get('AliJJetJtTask/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iFinder))
  hTrackJtTrue2D = f.Get('AliJJetJtTask/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iMCFinder))
  hTrackJtMisses2D = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackJtMisses2D/TrackJtMisses2DNFin{:02d}'.format(iFinder))
  hTrackJtFakes2D = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtUnfBg2D/JetConeJtUnfBg2DNFin{:02d}'.format(iFinder))
  
  hTrackMatchSuccess = [f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackMatchSuccess/TrackMatchSuccessNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)]
  
#   for h in hTrackJtCorr2D:
#     print("Bin 25702 content : {}".format(h.GetBinContent(25702)))
#     ybins = [h.GetYaxis().GetBinLowEdge(iBin) for iBin in range(1,h.GetNbinsY()+2)]
#     print(ybins)
#     ix = c_int()
#     iy = c_int()
#     iz = c_int()
#     h.GetBinXYZ(25702,ix,iy,iz)
#     print("{},{},{}".format(ix,iy,iz))
#     print("{}, Entries: {}".format(h.GetName(),h.GetEntries()))
#     for ibx in range(0,h.GetNbinsX()+1):
#       jtobs = h.GetXaxis().GetBinCenter(ibx)
#       for iby in range(0,h.GetNbinsY()+1):
#         ptobs = h.GetYaxis().GetBinCenter(iby)
#         for ibz in range(0,h.GetNbinsZ()+1):
#           jttrue = h.GetZaxis().GetBinCenter(ibz)
#           ib = h.GetBin(ibx,iby,ibz)
#           content = h.GetBinContent(ib)
#           #print("Bin {},{},{} is number {} (jtobs = {}, ptobs = {}, jttrue = {}, content = {})".format(ibx,iby,ibz,ib,jtobs,ptobs,jttrue,content))
#           if ib == 25702:
#             print("MOI")
#             return
#           if content > 0:
#             print("bin {},{},{} has content {}, (jtobs = {}, ptobs = {}, jttrue = {})".format(ibx,iby,ibz,content,jtobs,ptobs,jttrue))
#   return
  hTrackJtMeas2D.Add(hTrackJtFakes2D,-1)
  
  LogBinsJt = [hTrackJtMeas2D.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hTrackJtMeas2D.GetNbinsX()+2)]
  print(LogBinsJt)

  for h in hTrackJtCorr2D:
    print("{}".format(h.GetTitle()))
  hJetPtMeas = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iFinder))
  hJetPtMeasCoarse = Hist(jetBinBorders)
  for n,i in zip(numberJetsMeas,range(1,Njets+1)):
    hJetPtMeasCoarse.SetBinContent(i,n)
    
  hJetPtTrue = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iMCFinder))
  hJetPtTrueCoarse = Hist(jetBinBorders)
  for n,i in zip(numberJetsTrue,range(1,Njets+1)):
    hJetPtTrueCoarse.SetBinContent(i,n)
  LogBinsPt = [hJetPtTrue.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hJetPtTrue.GetNbinsX()+2)]
  hJetPtResponse = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorr/JetPtCorrNFin{:02d}'.format(iFinder))
  hJetPtResponseCoarse = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorrCoarse/JetPtCorrCoarseNFin{:02d}'.format(iFinder))


  TrackJtUnfolder = JtUnfolder.JtUnfolder('TrackJtUnfolder',jetBinBorders=jetBinBorders, Njets=Njets)
  TrackJtUnfolder.setTrackMatch(hTrackMatchSuccess)
  TrackJtUnfolder.drawTrackMatch("TrackMatch",'single')
  return
  TrackJtUnfolder.setPtBins(LogBinsPt)
  TrackJtUnfolder.setJtBins(LogBinsJt)
  TrackJtUnfolder.setJtMeas2D(hTrackJtMeas2D)
  TrackJtUnfolder.setJtTrue2D(hTrackJtTrue2D)
  TrackJtUnfolder.setMisses2D(hTrackJtMisses2D)
  TrackJtUnfolder.setFakes2D(hTrackJtFakes2D)
  TrackJtUnfolder.setJetPtMeas(hJetPtMeas)
  TrackJtUnfolder.setJetPtTrue(hJetPtTrue)
  TrackJtUnfolder.setJetPtMeasCoarse(hJetPtMeasCoarse)
  TrackJtUnfolder.setJetPtTrueCoarse(hJetPtTrueCoarse)
  TrackJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas, hJetPtResponse))
  TrackJtUnfolder.setJetPtResponseCoarse(createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse))
  TrackJtUnfolder.setNumberJetsMeas(numberJetsMeas)
  TrackJtUnfolder.setNumberJetsTrue(numberJetsTrue)
  TrackJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeas))
  TrackJtUnfolder.set2Dresponse(hTrackJtCorr2D)
  TrackJtUnfolder.unfold()
  #TrackJtUnfolder.plotResponse()
  TrackJtUnfolder.plotJetPt()
  TrackJtUnfolder.plotJt("PythiaTest.pdf")

  
  
if __name__ == "__main__": main()
  