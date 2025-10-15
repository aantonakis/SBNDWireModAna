#ifndef ANGLES_H
#define ANGLES_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;

// Return the new direction in rotated coordinates
void get_dir(float &xz, float &yz, int tpc, int plane, float trk_dirx, float trk_diry, float trk_dirz) 

{
  // First induction
  if (plane == 0) {
    // east
    if (tpc == 0) {
      float yp = trk_diry*TMath::Cos(TMath::Pi()/3) + trk_dirz*TMath::Sin(TMath::Pi()/3);              
      float zp = trk_diry*TMath::Sin(TMath::Pi()/3) - trk_dirz*TMath::Cos(TMath::Pi()/3);              
      yz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
      xz = TMath::ATan(trk_dirx/zp) * 180. / TMath::Pi();
    }
    // west
    else {
      float yp = trk_diry*TMath::Cos(TMath::Pi()/3) - trk_dirz*TMath::Sin(TMath::Pi()/3);              
      float zp = trk_diry*TMath::Sin(TMath::Pi()/3) + trk_dirz*TMath::Cos(TMath::Pi()/3);             
      yz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
      xz = TMath::ATan(trk_dirx/zp) * 180. / TMath::Pi();
    }
  }
  // 2nd induction
  else if (plane == 1) {
    // east
    if (tpc == 0) {
      float yp = trk_diry*TMath::Cos(TMath::Pi()/3) - trk_dirz*TMath::Sin(TMath::Pi()/3);              
      float zp = trk_diry*TMath::Sin(TMath::Pi()/3) + trk_dirz*TMath::Cos(TMath::Pi()/3);             
      yz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
      xz = TMath::ATan(trk_dirx/zp) * 180. / TMath::Pi();
    }
    // west
    else {
      float yp = trk_diry*TMath::Cos(TMath::Pi()/3) + trk_dirz*TMath::Sin(TMath::Pi()/3);              
      float zp = trk_diry*TMath::Sin(TMath::Pi()/3) - trk_dirz*TMath::Cos(TMath::Pi()/3);              
      yz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
      xz = TMath::ATan(trk_dirx/zp) * 180. / TMath::Pi();
    }
  }
  // Collection
  else {
    yz = TMath::ATan(trk_diry / trk_dirz) * 180. / TMath::Pi();
    xz = TMath::ATan(trk_dirx / trk_dirz) * 180. / TMath::Pi();
  }
  
}


#endif
