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


void profile_smooth_1d_itm(const char* input_file, const char* output_file, int DIM) {

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
          
            TH1D* h_result_temp = (TH1D*)h_summary->Clone(Form("h%d_%d", idx, DIM));
            h_result_temp->Reset();

            for (int i = 1; i < h_result_temp->GetNbinsX()+1; ++i) {
		std::string proj_str = num_str + "/projections/";
		proj_str += Form("proj_%d", i-1);
		const char* cstr_temp = num_str.c_str(); 
                //TH1D* h_1d_temp = (TH1D*)f->Get(cstr_temp+Form("proj_%d", i-1));
                TH1D* h_1d_temp = (TH1D*)f->Get(cstr_temp);
		TCanvas* c = new TCanvas("c", "c", 700, 500);
		h_1d_temp->Draw();
		c->Draw();
		gPad->WaitPrimitive();
		h_1d_temp->Draw();
                Float_t itm = 0.; // iterative truncated mean (ITM)
                Float_t itm_unc = 0.; // uncertainty of ITM
                Double_t itm_result[2];
                iterative_truncated_mean_std_err(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
                itm = itm_result[0];
                itm_unc = itm_result[1];
		h_result_temp->SetBinContent(i, itm);
		h_result_temp->SetBinError(i, itm_unc);

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
