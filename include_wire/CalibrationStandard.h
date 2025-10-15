#ifndef CALIBRATION_STANDARD_H
#define CALIBRATION_STANDARD_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"
#include "SCECorr.h"
#include "YZCorr.h"

using ROOT::Math::XYZVector;

const float calib_const_data0 = 0.0223037; 
const float calib_const_data1 = 0.0219534; 
const float calib_const_data2 = 0.0215156;

const float calib_const_mc0 = 0.0203521;
const float calib_const_mc1 = 0.0202351;
const float calib_const_mc2 = 0.0200727; 


float my_calib_const_corr(bool isData, int plane);

double dqdx_scale_correction_angle(double theta){

  double p0 = -1.55222;
  double p1 = 0.0500101;
  double p2 = -0.000211176;

  double this_bias = (p0 + p1 * theta + p2 * theta * theta) / 100.; // -- % to number
  double this_correction = 1. / (1. - this_bias);

  return this_correction;
}


XYZVector apply_sce_std(SCECorr *sce_corr, float &corr, int plane, const XYZVector& sp_sce_uncorr, float dirx, float diry, float dirz) 

{

    //XYZVector sp_sce_uncorr(sp_x, sp_y, sp_z);
    XYZVector sp_sce_corr = sce_corr -> WireToTrajectoryPosition(sp_sce_uncorr);
    double pitch_sce_uncorr = sce_corr -> meas_pitch(sp_sce_uncorr.X(), sp_sce_uncorr.Y(), sp_sce_uncorr.Z(), dirx, diry, dirz, plane, false);
    double pitch_sce_corr = sce_corr -> meas_pitch(sp_sce_uncorr.X(), sp_sce_uncorr.Y(), sp_sce_uncorr.Z(), dirx, diry, dirz, plane, true);
    corr = pitch_sce_uncorr / pitch_sce_corr;
    
    return sp_sce_corr;

}


void initialize_yz(YZCorr *yz_corr, bool isData) 

{
  //TString datapath = getenv("SBND_YZCORR_PATH");
  //std::cout << "YZ Map Path " << datapath << std::endl;  
  // Note: May need to update the maps in the future 
  TString yz_corr_f = "yz_correction_map_data1e20.root";
  if(!isData) yz_corr_f = "yz_correction_map_mcp2025b5e18.root";
  yz_corr -> SetFileStr(yz_corr_f);

}

// Calibration Constant Correction
float my_calib_const_corr(bool isData, int plane) 
{

  if (isData) return 1.;

  if (plane == 0) {
    return calib_const_data0/calib_const_mc0;
  }
  else if (plane == 1) {
    return calib_const_data1/calib_const_mc1;
  }
  else if (plane == 2) {
    return calib_const_data2/calib_const_mc2;
  }
  else {
    std::cout << "Invalid use of Calib Constant Correction!!!" << std::endl;
    return 1.;
  }
  return 1.;
}

#endif
