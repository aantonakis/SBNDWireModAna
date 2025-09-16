

// YZ Nonuniformity Calibrations
#include <TROOT.h>
#include <TH2F.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;



/*

  Nx = N planes and Ny = NTPCs
  Function loads correction hists for each plane/TPC 

*/
template <int NX, int NY>
void load_yz(TH2F (&hists)[NX][NY], const char* yz_file) {

    //TH2F* CzyHist_sce[kNplanes][2];

    TFile* file_SCEYZ = TFile::Open(yz_file, "READ");
    printf("wiremod_ndhist: Loading YZ nonuniformity correction histograms"
                " from %s\n", yz_file);
    
    for (int l = 0; l < NX; l++){
      for (int k = 0; k < NY; k++) {
        hists[l][k] = (TH2F*)file_SCEYZ->Get(Form("CzyHist_%i_%i",l,k));
      }
    }
    file_SCEYZ->Close();
    
    std::cout << hists[0][0] << "\n";

}

template <int NX, int NY>
double yz_corr(TH2F (&CzyHist_sce)[NX][NY], XYZVector sp, int ip, int tpc) {


  const int nbinz=100, nbinx=40, nbiny=80;
  const float lowz=0, highz=500, lowx=-200, highx=200, lowy=-200, highy=200;

  int ibinx = floor((sp.X() - lowx) / (highx - lowx) * nbinx);
  if (ibinx < 0 || ibinx >= nbinx)return 1;
  int ibiny = floor((sp.Y() - lowy) / (highy - lowy) * nbiny);
  if (ibiny < 0 || ibiny >= nbiny)return 1;
  int ibinz = floor((sp.Z() - lowz) / (highz - lowz) * nbinz);
  if (ibinz < 0 || ibinz >= nbinz)return 1;

  // YZ non-uniformity
  double CF_zy = CzyHist_sce[ip][tpc]->GetBinContent(ibinz+1, ibiny+1);

  return CF_zy;

}



