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
const UInt_t kNdims = 6;

const Float_t kTrackCut = 60.; // cm

const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "integral"};
const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", 
    "ThetaXZ (deg)", "ThetaYZ (deg)", "Integral"};

// Binnning --> Go as fine as possible and coarse grain later if needed
//const Int_t kNbins[kNdims] = { 460, 460, 560, 360, 360, 2500};
//const Double_t kXmin[kNdims] = { -230, -230, -30, -180, -180, 0};
//const Double_t kXmax[kNdims] = { 230, 230, 530, 180, 180, 5000};

// Copy for tracks minus one dimension
const Int_t kNbinsT[kNdims-1] = { 400, 400, 500, 72, 72};
const Double_t kXminT[kNdims-1] = { -200, -200, 0, -180, -180};
const Double_t kXmaxT[kNdims-1] = { 200, 200, 500, 180, 180};

SCECorr *sce_corr_mc = new SCECorr(false);
SCECorr *sce_corr_data = new SCECorr(true);


void ndmap_tracks_tpc_grid(TString list_file, TString out_suffix,

    // calibration options
    bool apply_sce = false,
    bool isData = false

) {

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

    TH1::AddDirectory(0);
 
    // 1 hist per plane per TPC. We also keep track of the number of tracks in
    // each eventual projection bin using TH2Is
    //THnSparseD* h[kNplanes * kNTPCs];
    THnSparseD* hTrack[kNplanes * kNTPCs];
    THnSparseD* hTrackFlag[kNplanes * kNTPCs]; // reset for each track

    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        //h[i] = new THnSparseD(Form("hwidth%d", i), "", kNdims, kNbins, kXmin, kXmax);
        hTrack[i] = new THnSparseD(Form("htrack%d", i), "", kNdims-1, kNbinsT, kXminT, kXmaxT);
        hTrackFlag[i] = new THnSparseD(Form("htrack%d", i), "", kNdims-1, kNbinsT, kXminT, kXmaxT);
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

      // Reset N-dimensional Track Counter
      for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        hTrackFlag[i]->Reset();
      }

      for (UInt_t ip = 0; ip < kNplanes; ip++) {
        // calculate the plane dependent angles
        float trk_thxz = -175.;
	float trk_thyz = -175.;
	unsigned short east_tpc = 0;

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
	  //
	  get_dir(trk_thxz, trk_thyz, my.tpc[ip][i], ip, *my.trk_dirx, *my.trk_diry, *my.trk_dirz);

          nevts++;

          XYZVector sp(my.x[ip][i], my.y[ip][i], my.z[ip][i]);

	  // ----------------- CALIBRATION BLOCK ------------------------ //
                    
	  XYZVector sp_sce;
	  float sce_q_corr = 1.;
		   
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
                   
          // ----------------- END CALIBRATION BLOCK ------------------------ //
                    

          // For Track Counting
          Double_t valT[kNdims-1]  = {
            sp_sce.X(), sp_sce.Y(), sp_sce.Z(), trk_thxz, trk_thyz
          };

          // select by TPC
          unsigned hit_idx = ip + kNplanes * my.tpc[ip][i];

          //h[hit_idx]->Fill(val);
          if (hTrackFlag[hit_idx]->GetBinContent(hTrackFlag[hit_idx]->GetBin(valT)) == 0) {
            hTrackFlag[hit_idx]->Fill(valT);
            hTrack[hit_idx]->Fill(valT);
          }
        } // loop over hits
      } // loop over planes
    } // loop over events
        
    //delete hTrackFlag;
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
      delete hTrackFlag[i];
    }
    
    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    
    std::cout << "About to write histograms to the output file" << std::endl;

    TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
    TString output_file_name = output_rootfile_dir + "/output_ndmap_tracks_tpc_" + out_suffix + ".root";
    out_rootfile = new TFile(output_file_name, "RECREATE");
    out_rootfile -> cd();
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
	std::cout << "Writing histograms for plane " << i << std::endl;
        hTrack[i]->Write();
	
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



