

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

const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void plot_contours(const char* input, const char* output_file) {

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


            TH2D* h_w = h->Projection(2, 0);
	    //int max = h_w->GetMaximum();
	    //int low = max/10

	    //Double_t levels[6] = {10, 100, 1000, 10000, 100000, 1000000};  // for example, contour at z = 50
            //h_w->SetContour(6, levels);
	    h_w->SetStats(0);
	    h_w->GetXaxis()->SetRangeUser(-90, 90);
	    h_w->GetXaxis()->SetTitle("Reconstructed #theta_{xz} [degrees]");
	    h_w->GetYaxis()->SetTitle("Hit Width [ticks]");
	    h_w->GetYaxis()->SetRangeUser(0, 10);
	    h_w->SetTitle(Form("TPC %d, Plane %d", tpc, plane));

	    cw->cd();
	    cw->SetLogz();
	    h_w->Draw("CONT Z LIST");
	    cw->Update();
	    cw->Write();	    
	   
        }

    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


