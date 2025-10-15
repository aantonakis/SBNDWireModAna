

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


void make_spline_1d(const char* input_mc, const char* input_data, const char* output_file, const char* label="x") {

    TFile* f_mc = TFile::Open(input_mc, "READ");
    if (!f_mc || f_mc->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }
    TFile* f_data = TFile::Open(input_data, "READ");
    if (!f_data || f_data->IsZombie()) {
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
            
	    std::string q_str = "hQ"+std::to_string(idx); // convert int to string
	    std::string w_str = "hW"+std::to_string(idx); // convert int to string

            const char* cstr_q = q_str.c_str(); 
            const char* cstr_w = w_str.c_str();
	    
            TH1D* h_q_mc = (TH1D*)f_mc->Get(cstr_q); 
            TH1D* h_q_data = (TH1D*)f_data->Get(cstr_q); 

            TH1D* h_w_mc = (TH1D*)f_mc->Get(cstr_w); 
            TH1D* h_w_data = (TH1D*)f_data->Get(cstr_w); 
            
            if (!h_q_mc) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f_mc->Close();
                f_data->Close();
                return;
            }

            TH1D* h_q_ratio = (TH1D*)h_q_data->Clone(Form("h_q_ratio_%d", idx));
            TH1D* h_w_ratio = (TH1D*)h_w_data->Clone(Form("h_w_ratio_%d", idx));
	
	    h_q_ratio->Divide(h_q_mc);
	    h_w_ratio->Divide(h_w_mc);
	   
	    h_q_ratio->GetXaxis()->SetTitle(label);
	    h_w_ratio->GetXaxis()->SetTitle(label);
	    
	    h_q_ratio->GetYaxis()->SetTitle("Integral [ADC]");
	    h_w_ratio->GetYaxis()->SetTitle("Width [ticks]");


	    // Create arrays for the spline input
    	    const int nbins = h_q_ratio->GetNbinsX();
    	    std::vector<double> xq, yq;
    	    std::vector<double> xw, yw;

            for (int i = 1; i <= nbins; ++i) {
              double binCenter = h_q_ratio->GetBinCenter(i);
              double binContentQ = h_q_ratio->GetBinContent(i);
              double binContentW = h_w_ratio->GetBinContent(i);
              if (binContentQ > 0) { // ignore empty bins if desired
                xq.push_back(binCenter);
                yq.push_back(binContentQ);
              }
              if (binContentW > 0) { // ignore empty bins if desired
                xw.push_back(binCenter);
                yw.push_back(binContentW);
              }
            }

            // Build a TSpline3 (cubic spline)
            TSpline3 *spline_q = new TSpline3(Form("temp_q_%d", idx), &xq[0], &yq[0], xq.size());
            TSpline3 *spline_w = new TSpline3(Form("temp_w_%d", idx), &xw[0], &yw[0], xw.size());
	    spline_q->SetName(Form("spline_q_%d", idx));
	    spline_w->SetName(Form("spline_w_%d", idx));
            spline_q->SetLineColor(kRed);
            spline_w->SetLineColor(kRed);

            spline_q->SetLineWidth(2);
            spline_w->SetLineWidth(2);

	    //spline_q->SetDirectory(0);
	    //spline_w->SetDirectory(0);

            outfile->cd();
           
            h_q_ratio->Write();
	    spline_q->Write(spline_q->GetName());	   

	    TCanvas* cq = new TCanvas(Form("c_q_%d", idx), "", 700, 500);
	    h_q_ratio->Draw("HISTE");
	    spline_q->Draw("Same");
	    //cq->Update();
	    cq->Write();	    

	    h_w_ratio->Write();	
	    spline_w->Write(spline_w->GetName());	   
	    
	    TCanvas* cw = new TCanvas(Form("c_w_%d", idx), "", 700, 500);
	    h_w_ratio->Draw("HISTE");
	    spline_w->Draw("Same");
	    //cw->Update();
	    cw->Write();	    

            delete h_q_mc;
            delete h_q_data;
            delete h_w_data;
            delete h_w_mc;
            delete h_q_ratio;
            delete h_w_ratio;
	    delete spline_q;
	    delete spline_w;
           
        
        }
    }
 
    delete f_mc;
    delete f_data;
    outfile->Close();

    printf("Finished processing.\n");

}


