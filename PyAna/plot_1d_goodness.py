import ROOT
import numpy as np
import sys


fm = ROOT.TFile.Open(sys.argv[1], "READ")
fd = ROOT.TFile.Open(sys.argv[2], "READ")


fm.ls()


for tpc in range(2):
	c = ROOT.TCanvas("c"+str(tpc), "", 700, 500)
	for plane in range(3):
		idx = 3*tpc + plane
		hgm = fm.Get("hG"+str(idx))  
		hgd = fd.Get("hG"+str(idx))
		hgd.GetXaxis().SetTitle("Reconstructed #theta_{xz}")
		hgd.GetYaxis().SetTitle("Goodness")
		hgd.GetXaxis().SetRangeUser(-90, 90)
		hgd.GetYaxis().SetRangeUser(0, 3)
		 
		if plane == 0:
			hgm.SetLineColor(4)
			hgd.SetLineColor(4)
			hgm.SetLineStyle(7)
			c.cd()
			hgd.Draw("HISTE")
			hgm.Draw("HISTE Same")
			c.Update()

		if plane == 1:
			hgm.SetLineColor(2)
			hgd.SetLineColor(2)
			hgm.SetLineStyle(7)
			c.cd()
			hgd.Draw("HISTE Same")
			hgm.Draw("HISTE Same")
			c.Update()
		if plane == 2:
			hgm.SetLineColor(1)
			hgd.SetLineColor(1)
			hgm.SetLineStyle(7)
			c.cd()
			hgd.Draw("HISTE Same")
			hgm.Draw("HISTE Same")
			c.Update()
	c.Draw()
	ROOT.gPad.WaitPrimitive()
	x = input("hello:")
	print(x)
fm.Close()
fd.Close()




