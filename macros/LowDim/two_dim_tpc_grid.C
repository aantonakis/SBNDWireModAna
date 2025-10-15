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


// Calibration Headers from Sungbin
#include "mylib.h"
#include "BetheBloch.h"
#include "SCECorr.h"
#include "YZCorr.h"

// Custom Helper Code
//#include "../../include/CalibrationStandard.h"
#include "CalibrationStandard.h"
#include "CalibNTupleInfo.h"
#include "Angles.h"

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);

//const UInt_t kNplanes = 3;
//const UInt_t kNTPCs = 2;
const UInt_t kNdims = 8;
const UInt_t kNdimsP = 5;

const Float_t kTrackCut = 60.; // cm

//const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "integral"};
//const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", 
//    "ThetaXZ (deg)", "ThetaYZ (deg)", "Integral"};

// Binnning --> Go as fine as possible and coarse grain later if needed
// x, y, z, txz, tyz, integral, width, goodness

const Int_t kNbins[kNdims] =   {  20,  20, 30, 20,  20, 1000, 200, 100};
const Double_t kXmin[kNdims] = { -200, -200, 0,  -180, -180, 0,    0,   0};
const Double_t kXmax[kNdims] = {  200,  200, 500, 180,  180, 3000, 20,  100};

BetheBloch *muon_BB = new BetheBloch(13); // setup for muons
SCECorr *sce_corr_mc = new SCECorr(false);
SCECorr *sce_corr_data = new SCECorr(true);
YZCorr *yz_corr = new YZCorr();
double lifetime = 100.; // mc default


void two_dim_tpc_grid(TString list_file, TString out_suffix,

    // calibration options
    bool apply_sce = false,
    bool apply_yz = false,
    bool apply_elife = false,
    bool apply_recom = false,
    bool isData = false,
    int dimx = 1,
    int dimy = 2

) {

    
    const Int_t kNbinsP[kNdimsP] = { kNbins[dimx],  kNbins[dimy], 1000, 200, 50};
    const Double_t kXminP[kNdimsP] = { kXmin[dimx], kXmin[dimy],  0,  0,  0};
    const Double_t kXmaxP[kNdimsP] = { kXmax[dimx], kXmax[dimy],  3000,  20,  100};


    // File List Management
    TChain *fChain = new TChain("caloskim/TrackCaloSkim");
    TString input_file_dir = getenv("DATA_PATH");
    TString sample_list_dir = getenv("SAMPLE_PATH");
    TString sample_list_label = getenv("FILELIST_LABEL");
    
    TString fileListPath = sample_list_dir + "/" + list_file;
    cout << "Opening : " << fileListPath << endl;

    std::ifstream file(fileListPath.Data());  // Convert TString to const char*
    if (!file) {
      cout << "File does not exist: " << fileListPath << endl;
      cout << "Exiting [ndhist_charges_tpc_crossers_grid]" << endl;
      return;
    }

    AddFilesToChain(fileListPath, fChain);
    MyCalib my(fChain);

    // SCE Calibration Initialization
    if (apply_sce) {
      if (isData) {
	sce_corr_data -> ReadHistograms();
      }
      else {
	sce_corr_mc -> ReadHistograms();
      }
    }
   
    // YZ Calibration Initialization
    if (apply_yz) {
      initialize_yz(yz_corr, isData);
      yz_corr -> ReadHistograms();
    }

    TH1::AddDirectory(0);
 
    // 1 hist per plane per TPC. We also keep track of the number of tracks in
    // each eventual projection bin using TH2Is
    THnSparseD* h[kNplanes * kNTPCs];

    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
      h[i] = new THnSparseD(Form("h2D%d", i), "", kNdimsP, kNbinsP, kXminP, kXmaxP);
    }

    size_t nevts = 0;
    size_t track_counter = 0;


    int track_idx = 0;
    while (my.reader.Next()) {
      track_idx++;
      if (*my.selected < 1) continue;

      // For CRT T0 study or TPC T0 study
      if (*my.whicht0 != 0) continue;

      // skip short tracks
      size_t nhits = my.rr[2].GetSize();
      if (nhits == 0) {
        fprintf(stderr, "Warning: Selected track (idx=%d, selected=%d) with no hits? Run=%d, Subrun=%d, Evt=%d. Skipping!\n", track_idx, *my.selected, *my.run, *my.subrun, *my.evt);
        continue;
      }
      if (my.rr[2][nhits - 1] < kTrackCut) continue;

      track_counter++;
            
      ROOT::Math::XYZVector trk_dir(*my.trk_dirx, *my.trk_diry, *my.trk_dirz);

      for (UInt_t ip = 0; ip < kNplanes; ip++) {
        // calculate the plane dependent angles
        float trk_thxz = -180.;
	float trk_thyz = -180.;

        for (size_t i = 0; i < my.x[ip].GetSize(); i++) {
          // skip nans
          if (my.x[ip][i] != my.x[ip][i]) continue;
                    
          // skip not on track
          if (!my.ontraj[ip][i]) continue;

          // goodness cut
          if (my.goodness[ip][i] >= 100.) continue;

          // hit trains have widths in increments of exactly 0.5
          // skip hits from these
          if (is_int(my.width[ip][i] * 2)) continue;
          //std::cout << "Current TPC " << tpc[ip][i] << std::endl; 
	  
	  // Angle code goes here! TODO
	  get_dir(trk_thxz, trk_thyz, my.tpc[ip][i], ip, *my.trk_dirx, *my.trk_diry, *my.trk_dirz);
          nevts++;

          XYZVector sp(my.x[ip][i], my.y[ip][i], my.z[ip][i]);

	  // ----------------- CALIBRATION BLOCK ------------------------ //
                    
	  XYZVector sp_sce;
	  float sce_q_corr = 1.;
	  float yz_q_corr = 1.;
	  float elife_q_corr = 1.;
	  float recom_q_corr = 1.;
		   
          if (apply_sce) {
            if (isData) {
	      sp_sce = apply_sce_std(sce_corr_data, sce_q_corr, ip, sp, *my.trk_dirx, *my.trk_diry, *my.trk_dirz);
	    }
	    else {
	      sp_sce = apply_sce_std(sce_corr_mc, sce_q_corr, ip, sp, *my.trk_dirx, *my.trk_diry, *my.trk_dirz);

            }
	  }
          else {
            sp_sce = sp;
	  }
	  if (apply_yz) {
            // Should probably be careful about using this without SCE corrections
            yz_q_corr = yz_corr -> GetYZCorr(sp_sce, ip);

	  }
	  if (apply_elife) {
	    if (isData) { 
	      if (my.tpc[ip][i] == 0) {	
    	        elife_q_corr = Lifetime_Correction(sp_sce.X(), 35.);
	      }
	      else {
	        elife_q_corr = Lifetime_Correction(sp_sce.X(), 35.);
	      }
	    }
	    else {
	      elife_q_corr = Lifetime_Correction(sp_sce.X(), lifetime);
            } 
	  }
          if (apply_recom) {
            recom_q_corr = my_calib_const_corr(isData, ip);
	  }
	  float total_q_corr = sce_q_corr * yz_q_corr * elife_q_corr * recom_q_corr;
                   

          // ----------------- END CALIBRATION BLOCK ------------------------ //
                    

	  float dimx_val = 0;
	  float dimy_val = 0;
	  if (dimx == 0) dimx_val = sp_sce.X();
	  if (dimx == 1) dimx_val = sp_sce.Y();
	  if (dimx == 2) dimx_val = sp_sce.Z();
	  if (dimx == 3) dimx_val = trk_thxz;
	  if (dimx == 4) dimx_val = trk_thyz;
	  
	  if (dimy == 0) dimy_val = sp_sce.X();
	  if (dimy == 1) dimy_val = sp_sce.Y();
	  if (dimy == 2) dimy_val = sp_sce.Z();
	  if (dimy == 3) dimy_val = trk_thxz;
	  if (dimy == 4) dimy_val = trk_thyz;

          Double_t val[kNdimsP]  = {
            dimx_val, dimy_val, my.integral[ip][i]*total_q_corr, my.width[ip][i], my.goodness[ip][i]
          };

          // select by TPC
          unsigned hit_idx = ip + kNplanes * my.tpc[ip][i];

          h[hit_idx]->Fill(val);
        } // loop over hits
      } // loop over planes
    } // loop over events
        

    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    
    std::cout << "About to write histograms to the output file" << std::endl;

    TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
    TString output_file_name = output_rootfile_dir + "/output_two_dim_tpc_" + out_suffix + ".root";
    out_rootfile = new TFile(output_file_name, "RECREATE");
    out_rootfile -> cd();
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
	std::cout << "Writing histograms for plane " << i << std::endl;
        h[i]->Write();
    }
   
    out_rootfile -> Close();

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



