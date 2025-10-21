

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
#include "../../../include_wire/Corr.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void plot_hit_trains(const char* input_file, const char* output_file) {
    
    gROOT->SetBatch(kTRUE);

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    
    
    for (int tpc = 0; tpc < kNTPCs; tpc++) {
        for (int plane = 0; plane < kNplanes; plane++) {
		
            // Get the histogram
            int idx = 3 * tpc + plane;
	    std::cout << "Plotting PLane " << idx << std::endl;
            //std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            std::string num_str = "h1D"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h = (THnSparseD*)f->Get(cstr);
            if (!h) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f->Close();
                return;
            }
	    TCanvas* cw = new TCanvas(Form("cw_large_%d", idx), "", 700, 500);
	    TLegend* leg = new TLegend(0.6, 0.4, 0.9, 0.9);
	    leg->SetBorderSize(0);
	    
	    bool start = false;
	    int color = 1;
    	    // Get ITM for each slice in the specified dimension
            for (int i = 1; i < h->GetAxis(0)->GetNbins() + 1; i++) {

	      if (i%4 != 0) continue;

              h->GetAxis(0)->SetRange(i, i);
              TH1D* h_1d_w = h->Projection(2); // width 
	      
	      float angle = h->GetAxis(0)->GetBinCenter(i);

	      h_1d_w->GetXaxis()->SetTitle("Width [ticks]");
	      h_1d_w->GetYaxis()->SetTitle("Area Normalized");
	      h_1d_w->SetTitle("");
	      h_1d_w->SetStats(0);
	      h_1d_w->SetLineWidth(2);

	     h_1d_w->Scale(1.0/h_1d_w->Integral());

	      if (angle > 70) {
		if (!start) {
	          h_1d_w->SetTitle(Form("TPC %d, Plane %d", tpc, plane));
		  h_1d_w->SetLineColor(color);      
		  h_1d_w->Draw("HISTE");
		  leg->AddEntry(h_1d_w, Form("#theta_{xz} = %.1f", angle));
		  start = true;
		  color += 1;
		}
		else {
		  if (color == 5) color += 1;
		  h_1d_w->SetLineColor(color);      
		  h_1d_w->Draw("HISTE Same");
                  leg->AddEntry(h_1d_w, Form("#theta_{xz} = %.1f", angle));
		  color += 1;
		}	
	     }
              
            }
	    leg->Draw("Same");

            outfile->cd();
	    cw->Write();

            delete h;
           
        }
    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


