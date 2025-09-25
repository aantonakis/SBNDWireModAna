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

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);

const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;
const UInt_t kNdims = 6;

const Float_t kTrackCut = 60.; // cm

const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "integral"};
const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", 
    "ThetaXZ (deg)", "ThetaYZ (deg)", "Integral"};

// Binnning --> Go as fine as possible and coarse grain later if needed
const Int_t kNbins[kNdims] = { 460, 460, 560, 360, 360, 2500};
const Double_t kXmin[kNdims] = { -230, -230, -30, -180, -180, 0};
const Double_t kXmax[kNdims] = { 230, 230, 530, 180, 180, 5000};

// Copy for tracks minus one dimension
const Int_t kNbinsT[kNdims-1] = { 460, 460, 560, 360, 360};
const Double_t kXminT[kNdims-1] = { -230, -230, -30, -180, -180};
const Double_t kXmaxT[kNdims-1] = { 230, 230, 530, 180, 180};

BetheBloch *muon_BB = new BetheBloch(13); // setup for muons
SCECorr *sce_corr_mc = new SCECorr(false);
SCECorr *sce_corr_data = new SCECorr(true);
YZCorr *yz_corr = new YZCorr();
double lifetime = 100.; // mc default


void ndhist_charges_tpc_crossers_grid_ntrack(TString list_file, TString out_suffix,

    // calibration options
    bool apply_sce = false,
    bool apply_yz = false,
    bool apply_elife = false,
    bool apply_recom = false,
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

    TTreeReader reader(fChain);

    if(isData) lifetime = 35.;

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
    THnSparseD* hTrack[kNplanes * kNTPCs];
    THnSparseD* hTrackFlag[kNplanes * kNTPCs]; // reset for each track

    TH2I* hi[kNplanes * kNTPCs * kNdims];
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i] = new THnSparseD(Form("hwidth%d", i), "", kNdims, kNbins, kXmin, kXmax);
        hTrack[i] = new THnSparseD(Form("htrack%d", i), "", kNdims-1, kNbinsT, kXminT, kXmaxT);
        hTrackFlag[i] = new THnSparseD(Form("htrack%d", i), "", kNdims-1, kNbinsT, kXminT, kXmaxT);
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j] = new TH2I(Form("hntrk_%d_%s", i, kLabels[j].Data()), "",
                    kNbins[j], kXmin[j], kXmax[j],
                    kNbins[kNdims - 1], kXmin[kNdims - 1], kXmax[kNdims - 1]);
        }
    }

    size_t nevts = 0;
    size_t track_counter = 0;


        //TTreeReader reader("caloskim/TrackCaloSkim", f);
        TTreeReaderValue<int> selected(reader, "trk.selected");
        TTreeReaderValue<int> whicht0(reader, "trk.whicht0");
        
        TTreeReaderValue<int> run(reader, "meta.run");
        TTreeReaderValue<int> evt(reader, "meta.evt");
        TTreeReaderValue<int> subrun(reader, "meta.subrun");

        // track length cut: use residual range on collection only
        //TTreeReaderArray<float> rr2(reader, "trk.hits2.rr");
        //TTreeReaderArray<float> rr1(reader, "trk.hits1.rr");
        //TTreeReaderArray<float> rr0(reader, "trk.hits0.rr");

        TTreeReaderValue<float> trk_dirx(reader, "trk.dir.x");
        TTreeReaderValue<float> trk_diry(reader, "trk.dir.y");
        TTreeReaderValue<float> trk_dirz(reader, "trk.dir.z");

        // this is super ugly, but pointers are worse for accessing TTreeReader below
        TTreeReaderArray<unsigned short> tpc[kNplanes] = {
            { reader, "trk.hits0.h.tpc" }, { reader, "trk.hits1.h.tpc" }, { reader, "trk.hits2.h.tpc" },
        };
        TTreeReaderArray<float> goodness[kNplanes] = {
            { reader, "trk.hits0.h.goodness" }, { reader, "trk.hits1.h.goodness" }, { reader, "trk.hits2.h.goodness" },
        };
        TTreeReaderArray<float> x[kNplanes] = {
            { reader, "trk.hits0.h.sp.x" }, { reader, "trk.hits1.h.sp.x" }, { reader, "trk.hits2.h.sp.x" },
        };
        TTreeReaderArray<float> y[kNplanes] = {
            { reader, "trk.hits0.h.sp.y" }, { reader, "trk.hits1.h.sp.y" }, { reader, "trk.hits2.h.sp.y" },
        };
        TTreeReaderArray<float> z[kNplanes] = {
            { reader, "trk.hits0.h.sp.z" }, { reader, "trk.hits1.h.sp.z" }, { reader, "trk.hits2.h.sp.z" },
        };
        TTreeReaderArray<float> width[kNplanes] = {
            { reader, "trk.hits0.h.width" }, { reader, "trk.hits1.h.width" }, { reader, "trk.hits2.h.width" },
        };
        TTreeReaderArray<float> integral[kNplanes] = {
            { reader, "trk.hits0.h.integral" }, { reader, "trk.hits1.h.integral" }, { reader, "trk.hits2.h.integral" },
        };
        TTreeReaderArray<bool> ontraj[kNplanes] = {
            { reader, "trk.hits0.ontraj" }, { reader, "trk.hits1.ontraj" }, { reader, "trk.hits2.ontraj" },
        };
        TTreeReaderArray<float> dirx[kNplanes] = {
            { reader, "trk.hits0.dir.x" }, { reader, "trk.hits1.dir.x" }, { reader, "trk.hits2.dir.x" },
        };
        TTreeReaderArray<float> diry[kNplanes] = {
            { reader, "trk.hits0.dir.y" }, { reader, "trk.hits1.dir.y" }, { reader, "trk.hits2.dir.y" },
        };
        TTreeReaderArray<float> dirz[kNplanes] = {
            { reader, "trk.hits0.dir.z" }, { reader, "trk.hits1.dir.z" }, { reader, "trk.hits2.dir.z" },
        };
        TTreeReaderArray<float> dqdx[kNplanes] = {
            { reader, "trk.hits0.dqdx" }, { reader, "trk.hits1.dqdx" }, { reader, "trk.hits2.dqdx" },
        };
        
        TTreeReaderArray<float> rr[kNplanes] = {
            { reader, "trk.hits0.rr" }, { reader, "trk.hits1.rr" }, { reader, "trk.hits2.rr" },
        };

        // keep track of how many tracks go into the hit distribution
        // since, e.g., all hits from 1 track will get the same angle, the
        // errors on angle will go like sqrt(NTracks) instead of sqrt(NHits)
        std::vector<TH2I*> hflag(kNplanes * kNTPCs * kNdims, nullptr);
        for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
            for (unsigned j = 0; j < kNdims; j++) {
                hflag.at(i * kNdims + j) = (TH2I*)hi[i * kNdims + j]->Clone(Form("hflag_%d_%d", i, j));
            }
        }

        int track_idx = 0;
        while (reader.Next()) {
            track_idx++;
            if (*selected < 1) continue;

            // For CRT T0 study or TPC T0 study
            if (*whicht0 != 0) continue;

            // skip short tracks
            size_t nhits = rr[2].GetSize();
            if (nhits == 0) {
                fprintf(stderr, "Warning: Selected track (idx=%d, selected=%d) with no hits? Run=%d, Subrun=%d, Evt=%d. Skipping!\n", track_idx, *selected, *run, *subrun, *evt);
                continue;
            }
            //if (rr2[nhits - 1] < kTrackCut) continue;
            if (rr[2][nhits - 1] < kTrackCut) continue;

            track_counter++;
            
            // TODO --> Need to Change to better coordinates
            ROOT::Math::XYZVector trk_dir(*trk_dirx, *trk_diry, *trk_dirz);
            //float trk_thxz = trk_dir.Theta() * 180. / TMath::Pi();
            //float trk_thyz = trk_dir.Phi() * 180. / TMath::Pi();

            // Reset N-dimensional Track Counter
            for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
                hTrackFlag[i]->Reset();
            }

            // reset track counting flags
            for (unsigned i = 0; i < kNplanes * kNTPCs * kNdims; i++) {
                hflag[i]->Reset();
            }

            for (UInt_t ip = 0; ip < kNplanes; ip++) {
              // calculate the plane dependent angles
              float trk_thxz = -175.;
	      float trk_thyz = -175.;
	      unsigned short east_tpc = 0;

                for (size_t i = 0; i < x[ip].GetSize(); i++) {
                    // skip nans
                    if (x[ip][i] != x[ip][i]) continue;
                    
                    // skip not on track
                    if (!ontraj[ip][i]) continue;

                    // goodness cut
                    if (goodness[ip][i] >= 100.) continue;

                    // hit trains have widths in increments of exactly 0.5
                    // skip hits from these
                    if (is_int(width[ip][i] * 2)) continue;
             	    //std::cout << "Current TPC " << tpc[ip][i] << std::endl; 
                    if (ip == 0) {
	              // First induction
	              if (tpc[ip][i] == east_tpc) {
		        // East TPC
	                float yp = *trk_diry*TMath::Cos(TMath::Pi()/3) + *trk_dirz*TMath::Sin(TMath::Pi()/3);              
	                float zp = *trk_diry*TMath::Sin(TMath::Pi()/3) - *trk_dirz*TMath::Cos(TMath::Pi()/3);              
		        trk_thyz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
                        trk_thxz = TMath::ATan(*trk_dirx/zp) * 180. / TMath::Pi();
		      }
		      else {
		        // West TPC --> Swap the Signs?
	                float yp = *trk_diry*TMath::Cos(TMath::Pi()/3) - *trk_dirz*TMath::Sin(TMath::Pi()/3);              
	                float zp = *trk_diry*TMath::Sin(TMath::Pi()/3) + *trk_dirz*TMath::Cos(TMath::Pi()/3);              
		        trk_thyz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
                        trk_thxz = TMath::ATan(*trk_dirx/zp) * 180. / TMath::Pi();
		      }
                    }
                    else if (ip == 1) {
                      // 2nd induction
	              if (tpc[ip][i] == east_tpc) {
		        // East TPC
	                float yp = *trk_diry*TMath::Cos(TMath::Pi()/3) - *trk_dirz*TMath::Sin(TMath::Pi()/3);              
	                float zp = *trk_diry*TMath::Sin(TMath::Pi()/3) + *trk_dirz*TMath::Cos(TMath::Pi()/3);              
		        trk_thyz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
                        trk_thxz = TMath::ATan(*trk_dirx/zp) * 180. / TMath::Pi();
		      }
		      else {
		        // West TPC --> Swap Signs?
	                float yp = *trk_diry*TMath::Cos(TMath::Pi()/3) + *trk_dirz*TMath::Sin(TMath::Pi()/3);              
	                float zp = *trk_diry*TMath::Sin(TMath::Pi()/3) - *trk_dirz*TMath::Cos(TMath::Pi()/3);              
		        trk_thyz = TMath::ATan(yp/zp) * 180. / TMath::Pi();
                        trk_thxz = TMath::ATan(*trk_dirx/zp) * 180. / TMath::Pi();
		      }
                    }
                    else {
                      // Collection Plane
		      trk_thyz = TMath::ATan(*trk_diry / *trk_dirz) * 180. / TMath::Pi();
                      trk_thxz = TMath::ATan(*trk_dirx / *trk_dirz) * 180. / TMath::Pi();

                    }

                    nevts++;

                    // for hit train study
                    /*
                    if (width[ip][i] > 10.4 && width[ip][i] < 10.7) {
                        printf("Unusually wide hit (w=%.1f, t=%.1f)! Run=%d, Evt=%d, Subrun=%d wire=%d, plane=%d TPC=%d\n",
                                width[ip][i], time[ip][i], *run, *evt, *subrun, wire[ip][i], ip, tpc[ip][i]);
                    }
                    */

                    XYZVector sp(x[ip][i], y[ip][i], z[ip][i]);

		    // ----------------- CALIBRATION BLOCK ------------------------ //
                    
		    XYZVector sp_sce;
		    float sce_q_corr = 1.;
		    float yz_q_corr = 1.;
		    float elife_q_corr = 1.;
		    float recom_q_corr = 1.;
		   
		    if (apply_sce) {
		      if (isData) {
		        sp_sce = apply_sce_std(sce_corr_data, sce_q_corr, ip, sp, *trk_dirx, *trk_diry, *trk_dirz);
		      }
		      else {
		        sp_sce = apply_sce_std(sce_corr_mc, sce_q_corr, ip, sp, *trk_dirx, *trk_diry, *trk_dirz);

		      }
		      //std::cout << "SP before SCE: X " << sp.X()     << " Y " << sp.Y()     << " Z " << sp.Z()     << std::endl;
		      //std::cout << "SP after  SCE: X " << sp_sce.X() << " Y " << sp_sce.Y() << " Z " << sp_sce.Z() << std::endl;
		      //std::cout << "SCE Corr " << sce_q_corr << std::endl;
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
			if (tpc[ip][i] == 0) {	
    		          elife_q_corr = Lifetime_Correction(sp_sce.X(), 44.5);
			}
			else {
			  elife_q_corr = Lifetime_Correction(sp_sce.X(), 33.8);
			}
		      }
		      else {
			elife_q_corr = Lifetime_Correction(sp_sce.X(), lifetime);
		      } 
		    }
	            
		    float total_q_corr = sce_q_corr * yz_q_corr * elife_q_corr * recom_q_corr;
                   
		    // TODO --> Ask Sunbin about this code block?

		        // trk angle dqdx reco bias corr: NOT APPLIED
		        //     //double this_dqdx_bias_corr = dqdx_scale_correction_angle(theta_trk_x);
		        //         //this_dqdx_bias_corr = 1.; // == no correction
		        //             //double dqdx_bias_corr = dqdx_all_corr * this_dqdx_bias_corr;

 
		    // ----------------- CALIBRATION BLOCK ------------------------ //
                    

		    //double dqdx_hit = dqdx[ip][i];
                    //double rr_hit = rr[ip][i];


		    // TODO --> This does not correct the track direction (maybe use start and end positions since these are through-going)
                    Double_t val[kNdims]  = {
                        sp_sce.X(), sp_sce.Y(), sp_sce.Z(), trk_thxz, trk_thyz, integral[ip][i]*total_q_corr 
                    };

                    // For Track Counting
                     Double_t valT[kNdims-1]  = {
                        sp_sce.X(), sp_sce.Y(), sp_sce.Z(), trk_thxz, trk_thyz
                    };


                    // select by TPC
                    unsigned hit_idx = ip + kNplanes * tpc[ip][i];
                    h[hit_idx]->Fill(val);
                    if (hTrackFlag[hit_idx]->GetBinContent(hTrackFlag[hit_idx]->FindFixBin(valT)) == 0) {
                        hTrackFlag[hit_idx]->Fill(valT);
                        hTrack[hit_idx]->Fill(valT);
                    }
                    // count track up to once per bin
                    for (unsigned j = 0; j < kNdims; j++) {
                        unsigned val_idx = hit_idx * kNdims + j;
                        Int_t target_bin = hflag[val_idx]->FindFixBin(val[j], val[kNdims - 1]);
                        if (hflag[val_idx]->GetBinContent(target_bin) == 0) {
                            hi[val_idx]->Fill(val[j], val[kNdims - 1]);
                            hflag[val_idx]->Fill(val[j], val[kNdims - 1]);
                        }
                        assert(hflag[val_idx]->GetBinContent(target_bin) != 0);
                    }
                } // loop over hits
            } // loop over planes
        } // loop over events
        
        for (unsigned i = 0; i < kNplanes * kNTPCs * kNdims; i++) {
            delete hflag.at(i);
        }
        
        delete hTrackFlag;
        //for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        //    delete hTrackFlag[i];
        //}
    

    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    

    TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
    TString output_file_name = output_rootfile_dir + "/output_ndhist_charges_tpc_crossers_" + out_suffix + ".root";
    out_rootfile = new TFile(output_file_name, "RECREATE");
    out_rootfile -> cd();
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i]->Write();
        hTrack[i]->Write();
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j]->Write();
        }
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



