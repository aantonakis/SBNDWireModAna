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

#include "../../include_wire/Fitting.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void profile_smooth_1d_landau(const char* input_file, const char* output_file, int DIM) {
    // Disable the web GUI (use the legacy display)
    //gROOT->SetBatch(kFALSE);
    //gEnv->SetValue("WebGui.HttpServer", "no");

    TFile* f = TFile::Open(input_file, "READ");
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
            std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            std::string summary_str = num_str + "/summary"; // append summary
            const char* cstr = summary_str.c_str(); 
            TH1D* h_summary = (TH1D*)f->Get(cstr);
            if (!h_summary) {
                std::cerr << "Histogram not found: " << cstr << "/summary" << std::endl;
                f->Close();
                return;
            }
	    //TCanvas* c1 = new TCanvas("c", "c", 700, 500);
            //h_summary->Draw();
            //fit->Draw("lsame");
	    //c1->Update();
	    //c1->WaitPrimitive();
          
            TH1D* h_result_temp = (TH1D*)h_summary->Clone(Form("h%d_%d", idx, DIM));
            h_result_temp->Reset();

            for (int i = 1; i < h_result_temp->GetNbinsX()+1; ++i) {
		std::string proj_str = num_str + "/projections/";
		proj_str += "proj_"; proj_str += std::to_string(i);
		//proj_str += Form("proj_%d", i);
		const char* cstr_temp = proj_str.c_str(); 
                //TH1D* h_1d_temp = (TH1D*)f->Get(cstr_temp+Form("proj_%d", i-1));
                TH1D* h_1d_temp = (TH1D*)f->Get(cstr_temp);
		
		// -------------- LANDAU FIT ---------------------- //
		
   		std::cout << "Fitting...\n";
		
		//const Double_t kQuantiles[2] = { 0.05, 0.999 };
		const Double_t kQuantiles[2] = { 0.04, 0.9 };
                Double_t limits[2];
                h_1d_temp->GetQuantiles(2, limits, kQuantiles);
 
   		// Setting fit range and start values
   		double fr[2];
   		double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
   		fr[0]=limits[0];
  	 	fr[1]=limits[1];
 
   		pllo[0]=0.; pllo[1]=0.; pllo[2]=1.0; pllo[3]=0.;
   		plhi[0]=1000; plhi[1]=2000; plhi[2]=10000000.0; plhi[3]=1000.0;
   		sv[0]=h_1d_temp->GetRMS()/2; sv[1]=h_1d_temp->GetMean(); sv[2]=50000.0; sv[3]=h_1d_temp->GetRMS()/2;
 
   		double chisqr;
   		int ndf;
   		TF1 *fit = langaufit(h_1d_temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
		
		if ((chisqr/ndf) > 5) {
		  std:cout << "Re-evaluate for a bad fit ..." << std::endl;
		  // probably a narrow gaus 
		  sv[0]=h_1d_temp->GetRMS()/3; sv[1]=h_1d_temp->GetMean(); sv[2]=50000.0; sv[3]=0.1;
		  const Double_t kQuantiles[2] = { 0.03, 0.9 };
                  Double_t limits[2];
                  h_1d_temp->GetQuantiles(2, limits, kQuantiles);
   		  fr[0]=limits[0];
  	 	  fr[1]=limits[1];
   		  fit = langaufit(h_1d_temp,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
		}
		

		std::cout << "Fitting done\nPlotting results...\n";
 
   		// Global style settings
   		gStyle->SetOptStat(1111);
   		gStyle->SetOptFit(111);
   		gStyle->SetLabelSize(0.03,"x");
   		gStyle->SetLabelSize(0.03,"y");
 
		//TCanvas* c2 = new TCanvas("c", "c", 700, 500);
		//h_1d_temp->Draw();
		//fit->SetLineWidth(5);
   		//fit->Draw("lsame");
		//c2->Update();
		//c2->WaitPrimitive();

		h_result_temp->SetBinContent(i, fp[1]);
		h_result_temp->SetBinError(i, fp[0]);

            }
	    h_result_temp->Write();
            // Perform the ITM calculation
            //profile_2d_proj_std_err(h_result_temp, input_file, DIM, tpc, plane, Ndim);

            //TH1D* hh = (TH1D*)h_result_temp->Clone(Form("h%d_%d", idx, DIM));

            //outfile->cd();
            //hh->Write();

            delete h_result_temp;
            delete h_summary;
           
        }
    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}
