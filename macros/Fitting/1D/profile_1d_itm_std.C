

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


void profile_1d_itm_std(const char* input_file, const char* output_file) {
    
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
            //std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            std::string num_str = "h1D"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h = (THnSparseD*)f->Get(cstr);
            if (!h) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f->Close();
                return;
            }
            TH1D* h_result_temp = h->Projection(0);
            h_result_temp->Reset();

            TH1D* h_charge = (TH1D*)h_result_temp->Clone(Form("hQ%d", idx));
            TH1D* h_width = (TH1D*)h_result_temp->Clone(Form("hW%d", idx));
            TH1D* h_good = (TH1D*)h_result_temp->Clone(Form("hG%d", idx));


	    TDirectory* P = outfile->mkdir(Form("Fits%d", idx));
	    TCanvas* cq_summary = new TCanvas(Form("cq_summary_%d", idx), "", 1400, 1000);
	    TCanvas* cw_summary = new TCanvas(Form("cw_summary_%d", idx), "", 1400, 1000);
	    
	    int nplots = h_result_temp->GetNbinsX();

	    int cols = TMath::CeilNint(TMath::Sqrt(nplots));
            int rows = TMath::CeilNint((double)nplots / cols);	    
	
	    cq_summary->Divide(cols, rows);
	    cw_summary->Divide(cols, rows);
	    std::vector<TH1D*> hqs;
	    std::vector<TH1D*> hws;

    	    // Get ITM for each slice in the specified dimension
            for (int i = 1; i < h->GetAxis(0)->GetNbins() + 1; i++) {
              h->GetAxis(0)->SetRange(i, i);
              TH1D* h_1d_q = h->Projection(1); // charge
              TH1D* h_1d_w = h->Projection(2); // width 
              TH1D* h_1d_g = h->Projection(3); // goodness
		
              Double_t itm_resultQ[2];
              Double_t itm_resultW[2];
              Double_t itm_resultG[2];

              iterative_truncated_mean_std_err(h_1d_q, -2, 1.75, 1.0e-4, itm_resultQ);
              iterative_truncated_mean_std_err(h_1d_w, -2, 1.75, 1.0e-4, itm_resultW);
              iterative_truncated_mean_std_err(h_1d_g, -2, 1.75, 1.0e-4, itm_resultG);
              
	      h_charge->SetBinContent(i, itm_resultQ[0]);
              h_charge->SetBinError(i, itm_resultQ[1]);
	      //h_charge->SetBinContent(i, silly_corr(plane)*itm_resultQ[0]);
              //h_charge->SetBinError(i, silly_corr(plane)*itm_resultQ[1]);
              h_width->SetBinContent(i, itm_resultW[0]);
              h_width->SetBinError(i, itm_resultW[1]);
              h_good->SetBinContent(i, itm_resultG[0]);
              h_good->SetBinError(i, itm_resultG[1]);
	
	      P->cd();
	      float edge_low = h->GetAxis(0)->GetBinCenter(i) - (h->GetAxis(0)->GetBinWidth(i))/2.;
	      float edge_high = h->GetAxis(0)->GetBinCenter(i) + (h->GetAxis(0)->GetBinWidth(i))/2;

	      h_1d_q->GetXaxis()->SetTitle("Integral [ADC]");
	      h_1d_w->GetXaxis()->SetTitle("Width [ticks]");

	      h_1d_q->SetTitle(Form("Bin Range: [%f, %f]", edge_low, edge_high));
	      h_1d_w->SetTitle(Form("Bin Range: [%f, %f]", edge_low, edge_high));

	      TCanvas* cq = new TCanvas(Form("cq_%d_%d", idx, i), "", 700, 500);
	      h_1d_q->Draw("HISTE");
	      TLine* lq = new TLine(itm_resultQ[0], 0, itm_resultQ[0], h_1d_q->GetMaximum());
	      lq->SetLineColor(2);
	      lq->SetLineWidth(2);
	      lq->Draw("Same");
	      cq->Write();
              
	      TCanvas* cw = new TCanvas(Form("cw_%d_%d", idx, i), "", 700, 500);
	      h_1d_w->Draw("HISTE");
	      TLine* lw = new TLine(itm_resultW[0], 0, itm_resultW[0], h_1d_w->GetMaximum());
	      lw->SetLineColor(2);
	      lw->SetLineWidth(2);
	      lw->Draw("Same");
	      cw->Write();

	      hqs.push_back(h_1d_q);		
	      hws.push_back(h_1d_w);		
	      
	      cq_summary->cd(i);
	      h_1d_q->Draw("HISTE");
	      //cq_summary->Update();	     
	      lq->Draw("Same");
	      //cq_summary->Update();	     
	      
	      cw_summary->cd(i);
	      h_1d_w->Draw("HISTE");
	      //cw_summary->Update();	     
	      lw->Draw("Same");
	      //cw_summary->Update();	     

	      //h_1d_q->Delete();
              //h_1d_w->Delete();
              //h_1d_g->Delete();
            }

	    P->cd();
	    cq_summary->Update();	     
	    cw_summary->Update();
	     
	    cq_summary->Write();
	    cw_summary->Write();
	    
	    hqs.clear();
	    hws.clear();

            outfile->cd();
            h_charge->Write();
            h_width->Write();
            h_good->Write();

            delete h_charge;
            delete h_width;
            delete h_good;
            delete h;
           
        }
    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


