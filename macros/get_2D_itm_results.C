

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

#include "itm_fit.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void get_2D_itm_results(const char* input_file, const char* output_file, 
    int DIMX, int DIMY, int Nx, int Ny, int Ndim) {
        

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
            std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h = (THnSparseD*)f->Get(cstr);
            if (!h) {
                std::cerr << "Histogram not found: hwidth" << idx << std::endl;
                f->Close();
                return;
            }


            TH2D* h_result_temp = h->Projection(DIMY, DIMX);
            h_result_temp->Reset();
  
            
            // Perform the ITM calculation
            profile_3d_proj(h_result_temp, input_file, DIMX, DIMY, Nx, Ny, tpc, plane, Ndim);

            TH2D* hh = (TH2D*)h_result_temp->Clone(Form("h_%d_%d_%d", idx, DIMX, DIMY));

            outfile->cd();
            hh->Write();

            delete h_result_temp;
            delete hh;
            delete h;
           
        }
    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


