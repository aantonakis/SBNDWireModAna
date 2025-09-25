/*
 * Create THnSparse containing hit widths and integrals along different dimensions
 * Dims: x, y, z, ThetaXZ, ThetaYZ
 * Input: Calibration ntuples. Select T0-tagged tracks from TPC
 * Gaus Hits: both widths and integrals
 * Option to apply different calibrations
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "Math/Vector3D.h"

#include "../../include_wire/Profile.h"

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);

const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;

void test_profile_1d(const char* input_file, const char* output_file,

    // minimum tracks per bin
    int minTracks = 100,
    int profileDim = 0,
    int projDim = 6

) 

{
   
    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }
    

    //TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    std::map<std::string, ProfileResult1D> results;

    // Get ND hists from the input file
    THnSparseD* h[kNplanes * kNTPCs];
    TH2I* hi[kNplanes * kNTPCs * kNdims];
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i] = f->Get(Form("hwidth%d", i));
       
        results[Form("hwidth%d", i)] = ProfileResult1D result = ProfileSparseDynamic1D(h[i],
          									       profileDim,
          									       projDim,
         									       minTracks);

    } 

    SaveAllProfileResult1D(results, output_file);
    delete f;    
}


std::vector<TString> filenames_from_input(const TString& input_arg, int nmax=-1) {
    std::vector<TString> filenames;
    if (input_arg.EndsWith(".root")) {
        filenames.push_back(input_arg);
        return filenames;
    }

    // read from file
    size_t nfiles = 0;
    std::ifstream ifile(input_arg);
    std::string line;
    while (std::getline(ifile, line)) {
        nfiles++;
        fprintf(stdout, "Adding file %zu: %s...\n", nfiles, line.c_str());
        filenames.push_back(TString(line));
        if (nfiles >= nmax && nmax > 0) break;
    }
    return filenames;
}


TString basename_prefix(const TString& input, const TString& prefix, const TString& suffix) {
    // remove path from filename and return new string with prefix or suffix added before extension
    TPRegexp re(".*/(.*)");
    TObjArray* matches = re.MatchS(input);
    TString result((static_cast<TObjString*>(matches->At(1)))->String());
    matches->Delete();
    return prefix + result;
}


bool is_int(Float_t val) {
    return std::abs(roundf(val) - val) < 0.00001f;
}



