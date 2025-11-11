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

#include "SelectionWire.h"

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);
bool is_one_third(float x, float tol = 1e-4);
bool is_two_thirds(float x, float tol = 1e-4);

const UInt_t kNdims = 10;

const Float_t kTrackCut = 60.; // cm

//const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "integral"};
const TString kTitles[kNdims] = { "x", "y", "z", "txz", "txy", "dqdx", "Q", "W", "G", "P"};
 
//    "ThetaXZ (deg)", "ThetaYZ (deg)", "Integral"};

// Binnning --> Go as fine as possible and coarse grain later if needed
// x, y, z, txz, tyz, integral, width, goodness

// Super Fine Width Binning for hit train testing
//
const Int_t kNbins[kNdims] =   { 200, 200, 250, 36, 36, 1000, 1000, 1600, 500, 2};

// The Limits are mostly finalized
const Double_t kXmin[kNdims] = { -200, -200, 0,  -90, -90, 0, 0,    0,   0,  0};
const Double_t kXmax[kNdims] = {  200,  200, 500, 90,  90, 3000, 3000, 16,  100, 2};

BetheBloch *muon_BB = new BetheBloch(13); // setup for muons
SCECorr *sce_corr_mc = new SCECorr(false);
SCECorr *sce_corr_data = new SCECorr(true);
YZCorr *yz_corr = new YZCorr();
double lifetime = 100.; // mc default

//const UInt_t kNdimsP = 4;

// Add 1 more for Pathological hits 
const UInt_t kNdimsP = 5;

const UInt_t kQ = 6; // Charge
const UInt_t kW = 7; // Width 
const UInt_t kG = 8; // Goodness

const UInt_t kP = 9; // Pathological Hits (small width at large angles)

void multi_dim_tracks_grid(TString list_file, TString out_suffix,

    // calibration options
    bool apply_sce = false,
    bool apply_yz = false,
    bool apply_elife = false,
    bool apply_recom = false,
    bool isData = false,
    std::vector<int> dim = {0},

    // Additional selections for special investigations 
    bool tpc_sel=true,
    bool crt_sel=false,
    bool pathological_sel=false,
    bool life_sel=false

) {

    std::cout << std::endl;    
    std::cout << "/---------------------------------------------------------------------------/" << std::endl;
    std::cout << std::endl;    
    std::cout << "Script Config:" << std::endl;
    std::cout << "Is this Data? " << isData << std::endl;
    std::cout << "Apply SCE: " << apply_sce << std::endl;
    std::cout << "Apply YZ: " << apply_yz << std::endl;
    std::cout << "Apply Lifetime: " << apply_elife << std::endl;
    std::cout << "Apply Calibration Const./Recomb.: " << apply_recom << std::endl;
    std::cout << std::endl;    
    std::cout << "Selections: " << std::endl;
    std::cout << "TPC Selection: " << tpc_sel << std::endl;
    std::cout << "CRT Selection: " << crt_sel << std::endl;
    std::cout << "Pathological Hit Selection: " << pathological_sel << std::endl;
    std::cout << "Lifetime Calibration Selection: " << life_sel << std::endl;
    std::cout << "DEBUG: kNplanes " << kNplanes << std::endl;
    std::cout << std::endl;    
    std::cout << "/---------------------------------------------------------------------------/" << std::endl;
    std::cout << std::endl;    


    // Add a pathological hit indicator at the end
    //const Int_t kNbinsP[kNdimsP] = { kNbins[dim],  kNbins[kQ], kNbins[kW], kNbins[kG], kNbins[kP]};
    //const Double_t kXminP[kNdimsP] = { kXmin[dim], kXmin[kQ], kXmin[kW], kXmin[kG], kXmin[kP]};
    //const Double_t kXmaxP[kNdimsP] = { kXmax[dim], kXmax[kQ], kXmax[kW], kXmax[kG], kXmax[kP]};


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
      cout << "Exiting [multi_dim_tpc_grid]" << endl;
      return;
    }

    AddFilesToChain(fileListPath, fChain);
    MyCalib my(fChain);

    // SCE Calibration Initialization
    if (apply_sce) {
      if (isData) {
	sce_corr_data -> ReadHistograms();
        // TODO DEBUG
        std::cout << "DATA DEBUG: Read SCE Hists" << std::endl;
      }
      else {
	sce_corr_mc -> ReadHistograms();
      }
    }
   
    // YZ Calibration Initialization
    if (apply_yz) {
      initialize_yz(yz_corr, isData);
      yz_corr -> ReadHistograms();
      std::cout << "DATA DEBUG: Initialized YZ" << std::endl;
    }

    TH1::AddDirectory(0);
 
    // 1 hist per plane per TPC. We also keep track of the number of tracks in
    // each eventual projection bin using TH2Is
    THnSparseD* h[kNplanes * kNTPCs];
    THnSparseD* hTracks[kNplanes * kNTPCs];
    THnSparseD* hTrackFlags[kNplanes * kNTPCs];

    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
      THnSparseD* h_temp = new THnSparseD(Form("h%d", i), "", kNdims, kNbins, kXmin, kXmax);
      THnSparseD* h_temp_trk = new THnSparseD(Form("hTrack%d", i), "", kNdims, kNbins, kXmin, kXmax);
      THnSparseD* h_temp_trk_flag = new THnSparseD(Form("hTrackFlags%d", i), "", kNdims, kNbins, kXmin, kXmax);
      //h[i] = new THnSparseD(Form("h1D%d", i), "", kNdimsP, kNbinsP, kXminP, kXmaxP);
      h[i] = static_cast<THnSparseD*>( h_temp->Projection(dim.size(), dim.data()) );
      hTracks[i] = static_cast<THnSparseD*>( h_temp_trk->Projection(dim.size(), dim.data()) );
      hTrackFlags[i] = static_cast<THnSparseD*>( h_temp_trk_flag->Projection(dim.size(), dim.data()) );
      h[i]->SetName(Form("hHit%d", i));
      hTracks[i]->SetName(Form("hTrack%d", i));

      // Set The axes labels       
      for (int j = 0; j < dim.size(); ++j) {
        h[i]->GetAxis(j)->SetTitle(kTitles[dim[j]]);
        hTracks[i]->GetAxis(j)->SetTitle(kTitles[dim[j]]);
      }

      h_temp->Delete();
      h_temp_trk->Delete();
      h_temp_trk_flag->Delete();
    }

    size_t nevts = 0;
    size_t track_counter = 0;

    int track_idx = 0;
    while (my.reader.Next()) {
      track_idx++;

      // Only use anode-cathode crossers for Lifetime study
      if ( (life_sel) && (*my.selected != 1) ) continue;

      // Main selections use both ACPTs and Cathode crossers
      if (*my.selected < 1) continue;

      // For TPC T0 study
      if ( (tpc_sel) && (*my.whicht0 != 0) ) continue;

      // CRT only T0 study      
      if ( (crt_sel) && (*my.whicht0 == 0) ) continue;
      
      // if neither TPC or CRT selection, then both are used 

      // skip short tracks
      size_t nhits = my.rr[2].GetSize();
      if (nhits == 0) {
        fprintf(stderr, "Warning: Selected track (idx=%d, selected=%d) with no hits? Run=%d, Subrun=%d, Evt=%d. Skipping!\n", track_idx, *my.selected, *my.run, *my.subrun, *my.evt);
        continue;
      }
      if (my.rr[2][nhits - 1] < kTrackCut) continue;

      track_counter++;

      // Reset N-dimensional Track Counter
      for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        hTrackFlags[i]->Reset();
      }
            
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
          
	  // This may kill the "pathological" hits
          if (my.mult[ip][i] > 1) continue;

          // hit trains have widths in increments of exactly 0.5
          // skip hits from these
          if (is_int(my.width[ip][i] * 2)) continue;
	  
	  // Angle code goes here!
	  get_dir(trk_thxz, trk_thyz, my.tpc[ip][i], ip, *my.trk_dirx, *my.trk_diry, *my.trk_dirz);
          
          // cut out large angles for lifetime correction
          if ( (life_sel) && (std::abs(trk_thxz) > 49) ) continue;

          // Hit trains also seem to come in thirds?	
          float frac = my.width[ip][i] - std::floor(my.width[ip][i]);
	  if (is_one_third(frac)) continue;
	  if (is_two_thirds(frac)) continue;

	  // ---------------------------------------------------------- //
	  //
	  // PATHOLOGICAL HIT Selection
	  //

          Double_t PATHOLOGICAL = 0.5; // --> Set to 0.5 for NOT pathological	  
          
	  unsigned IDX = ip + kNplanes * my.tpc[ip][i];
          bool cut_pathological = false;
	  try {
	    cut_pathological = txz_cut(trk_thxz, my.width[ip][i], IDX, isData); 
          }
          catch (const std::exception &e) {
            // Code to handle the error
            std::cerr << "Cut Pathological Error: " << e.what() << std::endl;
          }

	  if (cut_pathological) PATHOLOGICAL = 1.5;
	
	  // Can remove pathological hits if needed
          if ( (pathological_sel) && (cut_pathological) ) continue;	  
	  // ---------------------------------------------------------- //
          
          // TODO DATA DEBUG
          //if (isData) std::cout << "DATA DEBUG: Selected track --> entering calibration block" << std::endl;

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
          std::vector<double> vals;
          vals.reserve(dim.size());

          float dqdx_hit = my.dqdx[ip][i]*total_q_corr;

          for (int v = 0; v < dim.size(); ++v) {
	    double dim_val = 0;
	    if (dim[v] == 0) dim_val = sp_sce.X();
	    if (dim[v] == 1) dim_val = sp_sce.Y();
	    if (dim[v] == 2) dim_val = sp_sce.Z();
	    if (dim[v] == 3) dim_val = trk_thxz;
	    if (dim[v] == 4) dim_val = trk_thyz;
	    if (dim[v] == 5) dim_val = dqdx_hit;
	    if (dim[v] == 6) dim_val = my.integral[ip][i]*total_q_corr;
	    if (dim[v] == 7) dim_val = my.width[ip][i];
	    if (dim[v] == 8) dim_val = my.goodness[ip][i];
	    if (dim[v] == 9) dim_val = PATHOLOGICAL;
            vals.push_back(dim_val);
	  }

          // select by TPC
          unsigned hit_idx = ip + kNplanes * my.tpc[ip][i];

          // Fill the results
          h[hit_idx]->Fill(vals.data());
          if (hTrackFlags[hit_idx]->GetBinContent(hTrackFlags[hit_idx]->GetBin(vals.data())) == 0) {
            hTrackFlags[hit_idx]->Fill(vals.data());
            hTracks[hit_idx]->Fill(vals.data());
          }

        } // loop over hits
      } // loop over planes
    } // loop over events
   
    //delete hTrackFlag;
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
      delete hTrackFlags[i];
    }
    std::cout << "Finished the event loop ..." << std::endl;       

    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    
    std::cout << "About to write histograms to the output file" << std::endl;

    TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
    TString output_file_name = output_rootfile_dir + "/output_multi_dim_tracks_" + out_suffix + ".root";
    out_rootfile = new TFile(output_file_name, "RECREATE");
    out_rootfile -> cd();
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
	std::cout << "Writing histograms for plane " << i << std::endl;
        h[i]->Write();
        hTracks[i]->Write();
    }
   
    out_rootfile->Close();

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

bool is_one_third(float x, float tol = 1e-4) {
    return std::fabs(x - 1.0f/3.0f) < tol;
}

bool is_two_thirds(float x, float tol = 1e-4) {
    return std::fabs(x - 2.0f/3.0f) < tol;
}


