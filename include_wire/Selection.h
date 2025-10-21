#ifndef SELECTION_H
#define SELECTION_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;


// IMPORTANT: TRUE = Reject, FALSE = Keep
bool txz_cut(float txz, float width, int idx, bool isData){

  if (isData) {
    TFile* f = TFile::Open("CONST/data_txz_exclude.root", "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return true;
    }
    TF1* func = (TF1*)f->Get(Form("exclude_%d", idx));
    float y = func->Eval(txz);
    if (width < y) return true;
    return false;

  }
  else {
    TFile* f = TFile::Open("CONST/mc_txz_exclude.root", "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return true;
    }
    TF1* func = (TF1*)f->Get(Form("exclude_%d", idx));
    float y = func->Eval(txz);
    if (width < y) return true;
    return false;

  }

}

#endif
