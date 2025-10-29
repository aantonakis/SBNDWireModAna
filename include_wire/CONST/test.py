import ROOT
import sys

f = ROOT.TFile.Open(sys.argv[1], "READ")
f.ls()

f.Close()

