#ifndef SELECTION_WIRE_H
#define SELECTION_WIRE_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;


const double mc_fit[6][4] = {
    {2.0, 0.00022515895047805286, 2.0723916400591005e-08, 1.6688971261597966e-11},
    {2.0, 0.0002228709458460855, 2.713453903427959e-09, 1.9827225622990734e-11},
    {2.0, 0.00023251470884756188, -8.606819336303074e-09, 1.1971195470669035e-11},
    {2.0, 0.00022521659677347122, 2.030093110133196e-08, 1.69908580840407e-11},
    {2.0, 0.0002196124756275418, 7.4195000071279236e-09, 1.8543018840514705e-11},
    {2.0, 0.000233831543593686, -9.633043827975965e-09, 1.2107225865290702e-11}
};

const double data_fit[6][4] = {

    {2.0, 0.0001892230421399483, 1.4710503476323607e-08, 1.6064816317387788e-11},
    {2.0, 0.00017897498901686681, -6.492944516194252e-09, 1.7306266172753536e-11},
    {2.0, 0.0001346297699449016, 5.3104087995690904e-08, 2.0336015614461746e-12},
    {2.0, 0.00020000491090600708, 2.699342528050437e-10, 1.9085723056853728e-11},
    {2.0, 0.00019524018982115983, -6.413359042285074e-09, 1.7073518257560365e-11},
    {2.0, 0.00015509413830139263, 4.500935862965873e-08, 3.597088664976742e-12}
};

bool txz_cut(float txz, float width, unsigned idx, bool isData) {

  TF1 *func = new TF1(Form("fit_%d", idx), "[0] + [1]*x*x + [2]*x*x*x*x + [3]*x*x*x*x*x*x", -90, 90);
  if (isData) {
    for (int i = 0; i < 4; ++i) {
      func->SetParameter(i, data_fit[idx][i]);
    }
    float y = func->Eval(txz);
    delete func;
    if (width < y) return true;
    return false;
  }
  else {
    for (int i = 0; i < 4; ++i) {
      func->SetParameter(i, mc_fit[idx][i]);
    }
    float y = func->Eval(txz);
    delete func;
    if (width < y) return true;
    return false;
  }

}


// IMPORTANT: TRUE = Reject, FALSE = Keep

/*
bool txz_cut(float txz, float width, unsigned idx, bool isData){

  TString datapath = getenv("WIREMOD_WORKING_DIR");

  if (isData) {
    try {
      TFile* f = TFile::Open(datapath+"/include_wire/CONST/data_txz_exclude.root", "READ");
      if (!f || f->IsZombie()) {
        std::cerr << "File is a zombie mate" << std::endl;
        return true;
      }
      TF1* func = (TF1*)f->Get(Form("exclude_%d", idx));
      if (!func) {
        std::cerr << "Exclusion func didn't work!" << std::endl;
        return true;
      }
      float y = func->Eval(txz);
      if (width < y) return true;
      return false;
    }
    catch (const std::exception &e) {
      std::cerr << "Error in TXZ Cut: " << e.what() << std::endl;
      return false;
    }
  }
  else {
    try {
      TFile* f = TFile::Open(datapath+"/include_wire/CONST/mc_txz_exclude.root", "READ");
      if (!f || f->IsZombie()) {
        std::cerr << "File is a zombie mate" << std::endl;
        return true;
      }
      TF1* func = (TF1*)f->Get(Form("exclude_%d", idx));
      if (!func) {
        std::cerr << "Exclusion func didn't work!" << std::endl;
        return true;
      }
      float y = func->Eval(txz);
      if (width < y) return true;
      return false;
    }   
    catch (const std::exception &e) {
      std::cerr << "Error in TXZ Cut: " << e.what() << std::endl;
      return false;
    }
  }

}
*/


#endif
