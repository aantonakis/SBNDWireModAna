#ifndef CORR_H
#define CORR_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;

const float calib_const_data0 = 0.0223037; 
const float calib_const_data1 = 0.0219534; 
const float calib_const_data2 = 0.0215156;

const float calib_const_mc0 = 0.0203521;
const float calib_const_mc1 = 0.0202351;
const float calib_const_mc2 = 0.0200727; 


//float silly_corr(bool isData, int plane);

// Calibration Constant Correction
float silly_corr(int plane) 
{

  //if (isData) return 1.;

  if (plane == 0) {
    return (calib_const_data0/calib_const_mc0)*(calib_const_data0/calib_const_mc0);
  }
  else if (plane == 1) {
    return (calib_const_data1/calib_const_mc1)*(calib_const_data1/calib_const_mc1);
  }
  else if (plane == 2) {
    return (calib_const_data2/calib_const_mc2)*(calib_const_data2/calib_const_mc2);
  }
  else {
    std::cout << "Invalid use of Calib Constant Correction!!!" << std::endl;
    return 1.;
  }
  return 1.;
}

#endif
