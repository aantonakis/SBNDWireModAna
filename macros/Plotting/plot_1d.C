

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


void plot_1d(const char* input, const char* output_file, int dim=0) {

    TFile* f = TFile::Open(input, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);
    
    for (int tpc = 0; tpc < kNTPCs; tpc++) {
	
	TCanvas* cq = new TCanvas(Form("c_q_%d", tpc), "", 700, 500);
	TCanvas* cw = new TCanvas(Form("c_w_%d", tpc), "", 700, 500);

        for (int plane = 0; plane < kNplanes; plane++) {

            // Get the histogram
            int idx = 3 * tpc + plane;

            TH1D* h_q_ratio = (TH1D*)f->Get(Form("h_q_ratio_%d", idx));
            TH1D* h_w_ratio = (TH1D*)f->Get(Form("h_w_ratio_%d", idx));
            
	    TSpline3* spline_q = (TSpline3*)f->Get(Form("spline_q_%d", idx));
	    TSpline3* spline_w = (TSpline3*)f->Get(Form("spline_w_%d", idx));
	    spline_q->SetLineStyle(7);
	    spline_w->SetLineStyle(7);

	    h_q_ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	    h_w_ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	    
	    h_q_ratio->GetYaxis()->SetTitle("Hit Charge: Data/Simulation");
	    h_w_ratio->GetYaxis()->SetTitle("Hit Width: Data/Simulation");

	    if (dim == 0) {
	     if (tpc == 0) {
	       h_q_ratio->GetXaxis()->SetRangeUser(-200, 0);
	       h_w_ratio->GetXaxis()->SetRangeUser(-200, 0);
	     }
	     else {
	       h_q_ratio->GetXaxis()->SetRangeUser(0, 200);
	       h_w_ratio->GetXaxis()->SetRangeUser(0, 200);
	     }
	    } 
	    if (dim == 3) {
	       h_q_ratio->GetXaxis()->SetRangeUser(-90, 90);
	       h_w_ratio->GetXaxis()->SetRangeUser(-90, 90);
	    }

	    if (plane == 0) {
	      h_q_ratio->SetLineColor(4);
	      h_w_ratio->SetLineColor(4);
	      spline_q->SetLineColor(kBlue);
	      spline_w->SetLineColor(kBlue);
	
	      h_q_ratio->SetTitle(Form("TPC %d", tpc));
	      h_q_ratio->SetStats(0);

	      h_w_ratio->SetTitle(Form("TPC %d", tpc));
	      h_w_ratio->SetStats(0);

	      cq->cd();
	      h_q_ratio->Draw("HISTE");
	      spline_q->Draw("Same");
	      cq->Update();
		

	      cw->cd();
	      h_w_ratio->Draw("HISTE");
	      spline_w->Draw("Same");
	      cw->Update();
	      

	      //TLegend* leg0 = new TLegend(0.5, 0.1, 0.);

	    }
	    else if (plane == 1) {
	      h_q_ratio->SetLineColor(2);
	      h_w_ratio->SetLineColor(2);

	      cq->cd();
	      h_q_ratio->Draw("HISTE Same");
	      spline_q->Draw("Same");
	      cq->Update();
	
	      cw->cd();
	      h_w_ratio->Draw("HISTE Same");
	      spline_w->Draw("Same");
	      cw->Update();

	    }
	    else {
	      h_q_ratio->SetLineColor(1);
	      h_w_ratio->SetLineColor(1);
	      spline_q->SetLineColor(kBlack);
	      spline_w->SetLineColor(kBlack);

	      cq->cd();
	      h_q_ratio->Draw("HISTE Same");
	      spline_q->Draw("Same");
	      cq->Update();
	
	      cw->cd();
	      h_w_ratio->Draw("HISTE Same");
	      spline_w->Draw("Same");
	      cw->Update();
	      
	    }

            //delete h_q_ratio;
            //delete h_w_ratio;
	    //delete spline_q;
	    //delete spline_w;
           
        
        }
	//cq->Update();
	//cw->Update();

        outfile->cd();   

	cq->Write();	    
	cw->Write();	    

    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


