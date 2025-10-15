

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "Math/Vector3D.h"

#include "../../../include_wire/Fitting.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void profile_1d_itm_std_flist(const char* filelist, const char* output_file) {

    std::ifstream infile(filelist);
    if (!infile.is_open()) {
        Error("SumHistograms", "Cannot open file list: %s", filelist);
        return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    // initialize the array of multidimensional hists
    THnSparseD* h[kNplanes * kNTPCs] = {nullptr};
 
    std::string fname;
    int nfiles = 0;

    while (infile >> fname) {
        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            Warning("SumHistograms", "Skipping bad file: %s", fname.c_str());
            continue;
        }
	// Get all the hists for this file 
        for (int tpc = 0; tpc < kNTPCs; tpc++) {
          for (int plane = 0; plane < kNplanes; plane++) {
            int idx = 3 * tpc + plane;
            
            std::string num_str = "h1D"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h_temp = (THnSparseD*)f->Get(cstr);
            if (!h_temp) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f->Close();
                continue;
            }
	    // Initialize accumulator for this index if needed
            if (!h[idx]) {
                h[idx] = dynamic_cast<THnSparseD*>(h_temp->Clone(Form("h1DNew%d", idx)));
                //h[idx]->SetDirectory(nullptr);  // Detach from input file
                h[idx]->Reset();                // Zero contents, keep binning
            }

	    h[idx]->Add(h_temp);
          }
        }
        nfiles++;
        f->Close();
    }


    for (int tpc = 0; tpc < kNTPCs; tpc++) {
        for (int plane = 0; plane < kNplanes; plane++) {

            // Get the histogram
            int idx = 3 * tpc + plane;
            TH1D* h_result_temp = h[idx]->Projection(0);
            h_result_temp->Reset();

            TH1D* h_charge = (TH1D*)h_result_temp->Clone(Form("hQ%d", idx));
            TH1D* h_width = (TH1D*)h_result_temp->Clone(Form("hW%d", idx));
            TH1D* h_good = (TH1D*)h_result_temp->Clone(Form("hG%d", idx));

    	    // Get ITM for each slice in the specified dimension
            for (int i = 1; i < h[idx]->GetAxis(0)->GetNbins() + 1; i++) {
              h[idx]->GetAxis(0)->SetRange(i, i);
              TH1D* h_1d_q = h[idx]->Projection(1); // charge
              TH1D* h_1d_w = h[idx]->Projection(2); // width 
              TH1D* h_1d_g = h[idx]->Projection(3); // goodness
		
              Double_t itm_resultQ[2];
              Double_t itm_resultW[2];
              Double_t itm_resultG[2];

              iterative_truncated_mean_std_err(h_1d_q, -2, 1.75, 1.0e-4, itm_resultQ);
              iterative_truncated_mean_std_err(h_1d_w, -2, 1.75, 1.0e-4, itm_resultW);
              iterative_truncated_mean_std_err(h_1d_g, -2, 1.75, 1.0e-4, itm_resultG);
              
	      h_charge->SetBinContent(i, itm_resultQ[0]);
              h_charge->SetBinError(i, itm_resultQ[1]);
              h_width->SetBinContent(i, itm_resultW[0]);
              h_width->SetBinError(i, itm_resultW[1]);
              h_good->SetBinContent(i, itm_resultG[0]);
              h_good->SetBinError(i, itm_resultG[1]);
              
	      h_1d_q->Delete();
              h_1d_w->Delete();
              h_1d_g->Delete();
            }

            outfile->cd();
            h_charge->Write();
            h_width->Write();
            h_good->Write();

            delete h_charge;
            delete h_width;
            delete h_good;
            //delete h;
           
        }
    }
    for (unsigned i = 0; i < 3 * kNTPCs; i++) {
      delete h[i];
    }

    //delete h;
    //delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


