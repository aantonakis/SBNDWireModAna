/*
 * Create THnSparse containing hit widths along different dimensions
 * x, y, z, ThetaXZ, ThetaYZ, dQdx, residual range
 * Input: Calibration ntuples. Select T0-tagged tracks
 * No Calibration in this Module !!!
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

//#include "SCECorr.h"

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);

const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;
const UInt_t kNdims = 7;
const Float_t kTrackCut = 60.; // cm

const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "dqdx", "integral" };
const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", 
    "ThetaXZ (deg)", "ThetaYZ (deg)", "dQ/dx", "Integral (ADC)" };
const Int_t kNbins[kNdims] = { 600, 600, 700, 360, 360, 50, 2500 };
const Double_t kXmin[kNdims] = { -300, -300, -100, -180, -180, 0, 0 };
const Double_t kXmax[kNdims] = { 300, 300, 600, 180, 180, 5000, 5000 };



void wiremod_ndintegrals_nocalib_fine(const char* input_file, const char* output_file) {

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);
 
    // 1 hist per plane per TPC. We also keep track of the number of tracks in
    // each eventual projection bin using TH2Is
    THnSparseD* h[kNplanes * kNTPCs];
    TH2I* hi[kNplanes * kNTPCs * kNdims];
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i] = new THnSparseD(Form("hwidth%d", i), "", kNdims, kNbins, kXmin, kXmax);
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j] = new TH2I(Form("hntrk_%d_%s", i, kLabels[j].Data()), "",
                    kNbins[j], kXmin[j], kXmax[j],
                    kNbins[kNdims - 1], kXmin[kNdims - 1], kXmax[kNdims - 1]);
        }
    }

    size_t nevts = 0;
    size_t track_counter = 0;


        TTreeReader reader("caloskim/TrackCaloSkim", f);
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
            // Not Stopping Muons
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
                    double dqdx_hit = dqdx[ip][i];
                    //double rr_hit = rr[ip][i];

                    Double_t val[kNdims]  = {
                        sp.X(), sp.Y(), sp.Z(), trk_thxz, trk_thyz, dqdx_hit, integral[ip][i] 
                    };

                    // select by TPC
                    unsigned hit_idx = ip + kNplanes * tpc[ip][i];
                    h[hit_idx]->Fill(val);

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

        delete f;
    

    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    
    outfile->cd();   
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i]->Write();
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j]->Write();
        }
    }

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



