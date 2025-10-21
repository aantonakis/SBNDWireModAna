

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


void plot_exclusion(const char* input, const char* output_file) {

    gROOT->SetBatch(kTRUE);
    
    TFile* f = TFile::Open(input, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);
    
    outfile->cd();   
    for (int tpc = 0; tpc < kNTPCs; tpc++) {
	
        for (int plane = 0; plane < kNplanes; plane++) {
	
            // Get the histogram
            int idx = 3 * tpc + plane;
	    TCanvas* cw = new TCanvas(Form("c_w_%d", idx), "", 700, 500);
            std::string num_str = "h1D"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h = (THnSparseD*)f->Get(cstr);
            if (!h) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f->Close();
                return;
            }

            TH1D* h_result_width = h->Projection(0);
            h_result_width->Reset();
            
	    TH2D* h_w = h->Projection(2, 0);
	    h_w->SetStats(0);
	    h_w->GetXaxis()->SetRangeUser(-90, 90);
	    h_w->GetXaxis()->SetTitle("Reconstructed #theta_{xz} [degrees]");
	    h_w->GetYaxis()->SetTitle("Hit Width [ticks]");
	    h_w->GetYaxis()->SetRangeUser(0, 10);
	    h_w->SetTitle(Form("TPC %d, Plane %d", tpc, plane));


	    // Loop over the bins to define the MPV curve
            for (int i = 1; i < h->GetAxis(0)->GetNbins() + 1; i++) {
	      float angle = h->GetAxis(0)->GetBinCenter(i);
	      if ( (angle > 60) || (angle < -60) ) continue;

              h->GetAxis(0)->SetRange(i, i);
	      
              TH1D* h_1d_w = h->Projection(2); // width 
		
              Double_t itm_resultW[2];

              iterative_truncated_mean_std_err(h_1d_w, -2, 1.75, 1.0e-4, itm_resultW);
              
	      h_result_width->SetBinContent(i, itm_resultW[0]);
	      h_result_width->SetBinError(i, itm_resultW[1]);

	    }

	    h_result_width->SetLineColor(2);
	    h_result_width->SetLineWidth(2);

	    cw->cd();
	    cw->SetLogz();
	    h_w->Draw("Colz");
	    h_result_width->Draw("HISTE Same");
	    cw->Update();
	    cw->Write();	    
	   
        }

    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


