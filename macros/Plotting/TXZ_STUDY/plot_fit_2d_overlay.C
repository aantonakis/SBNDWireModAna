

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


void plot_fit_2d_overlay(const char* input, const char* output_file, const char* xlabel = "Reconstructed #theta_{xz} [degrees]") {

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
	    h_w->GetXaxis()->SetTitle(xlabel);
	    h_w->GetYaxis()->SetTitle("Hit Width [ticks]");
	    h_w->GetYaxis()->SetRangeUser(0, 10);
	    h_w->SetTitle(Form("TPC %d, Plane %d", tpc, plane));


	    // Loop over the bins to define the MPV curve
            for (int i = 1; i < h->GetAxis(0)->GetNbins() + 1; i++) {

              h->GetAxis(0)->SetRange(i, i);
	      
              TH1D* h_1d_w = h->Projection(2); // width 
		
              Double_t itm_resultW[2];

              iterative_truncated_mean_std_err(h_1d_w, -2, 1.75, 1.0e-4, itm_resultW);
              
	      h_result_width->SetBinContent(i, itm_resultW[0]);
	      h_result_width->SetBinError(i, itm_resultW[1]);

	    }

	    // Define a quadratic function: f(x) = p0 + p1*x + p2*x^2
            double xmin = -60.0;
            double xmax = 60.0;
            TF1 *func = new TF1(Form("fit_%d", idx), "[0] + [1]*x*x + [2]*x*x*x*x + [3]*x*x*x*x*x*x", xmin, xmax);

	    /*
	    TF1 *fexp = new TF1("fexp",
                        "[0]*exp([1]*abs(x)) + [2]",
                        xmin, xmax);
	    */

            // (Optional) Set initial parameter guesses
            func->SetParameters(3.5, 1e-4, 1e-8, 1e-8);

            // Fit histogram in the specified range
            h_result_width->Fit(func, "R");  // "R" = use range defined in TF1

	    func->SetRange(-90, 90);
	    TF1* func2 = (TF1*)func->Clone(Form("exclude_%d", idx));
	    func2->SetParameter(0, 2);
	    func2->SetLineColor(2);
	  	   

	    h_result_width->SetLineColor(2);
	    h_result_width->SetLineWidth(4);

	    cw->cd();
	    cw->SetLogz();
	    h_w->Draw("Colz");
	    h_result_width->Draw("HISTE Same");
	    func->SetLineColor(4);
	    func->Draw("Same");
	    func2->Draw("Same");
	    cw->Update();
	    cw->Write();	    
	    func2->Write();
	   
        }

    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


