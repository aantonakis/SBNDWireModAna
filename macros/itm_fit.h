// Header File for itm_fit.cc

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

void hist_mean_unc(const TH1* h, Double_t tol, Double_t* result);
void iterative_truncated_mean(TH1* h, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result);
void profile_2d_proj(TH1D* h_out, const char* input_file, int dim, int tpc, int plane, int Ndim);
void profile_3d_proj(TH2D* h_out, const char* input_file, int dimx, int dimy, int Nx, int Ny, int tpc, int plane, int Ndim);


/*

Iterative truncated mean calculation for histograms that uses bootstrapping to estimate the uncertainty 


*/


// return the mean & uncertainty on the mean based on bins in a histogram, accounting for each bin's error
void hist_mean_unc(const TH1* h, Double_t tol, Double_t* result) {
    Int_t nbins = h->GetNbinsX();
    if (h->Integral() == 0) return; 

    const Int_t kNthrows = 1000;
    Double_t means[kNthrows];

    TRandom3 rng;
    for (Int_t i = 0; i < kNthrows; i++) {
        TH1* htmp = (TH1*)h->Clone();
        htmp->Reset();
        for (Int_t k = 1; k <= nbins; k++) {
            Double_t bc = h->GetBinContent(k) + rng.Gaus(0, h->GetBinError(k));
            htmp->SetBinContent(k, TMath::Max(bc, 0.));
        }
        means[i] = htmp->GetMean();
        delete htmp;
    }

    result[0] = TMath::Mean(kNthrows, means);
    result[1] = TMath::RMS(kNthrows, means);
}


void iterative_truncated_mean(TH1* h, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result) {
    if (sig_down > sig_up) {
        fprintf(stderr, "Warning: reversing iterative truncated mean limits [%.2e,%.2e]\n", sig_down, sig_up);
        Double_t tmp = sig_down;
        sig_down = sig_up;
        sig_up = tmp;
    }

    Double_t mean = h->GetMean();
    Double_t sd = h->GetRMS();
    std::cout << mean << ", " << sd << ", " << result[0] << "\n";

    // return mean, std err of iterative truncated mean
    if (TMath::Abs(result[0] - mean) < tol) {
        hist_mean_unc(h, 1.0e-6, result);
        return;
    }
    result[0] = mean;
    result[1] = sd;
    
    // get positions to cut based on quantiles
    Double_t probs[1] = { 0.5 };
    Double_t median[1];
    h->GetQuantiles(1, median, probs);

    std::cout << "median: " << median[0] << "\n";
    std::cout << median[0] + sig_down * sd << ", " << median[0] + sig_up * sd << "\n";
    
    TH1* hnew = (TH1*)h->Clone();
    hnew->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinLowEdge(i) > median[0] + sig_up * sd || h->GetBinLowEdge(i) + h->GetBinWidth(i) < median[0] + sig_down * sd) continue;
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }
    
    iterative_truncated_mean(hnew, sig_down, sig_up, tol, result);
    delete hnew;
    return;
}


void profile_2d_proj(TH1D* h_out, const char* input_file, int dim, int tpc, int plane, int Ndim) {
    // Open the input file
    TFile* rfile = TFile::Open(input_file);
    if (!rfile || rfile->IsZombie()) {
        std::cerr << "Error opening file: " << input_file << std::endl;
        return;
    }

    // Get the histogram
    int idx = 3 * tpc + plane;
    std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
    const char* cstr = num_str.c_str();
    THnSparseD* h = (THnSparseD*)rfile->Get(cstr);
    if (!h) {
        std::cerr << "Histogram not found: hwidth" << idx << std::endl;
        rfile->Close();
        return;
    }

    // Get ITM for each slice in the specified dimension
    for (int i = 1; i < h->GetAxis(dim)->GetNbins() + 1; i++) {
        h->GetAxis(dim)->SetRange(i, i);
        TH1D* h_1d_temp = h->Projection(Ndim -1);
        Float_t itm = 0.; // iterative truncated mean (ITM)
        Float_t itm_unc = 0.; // uncertainty of ITM
        Double_t itm_result[2];
        iterative_truncated_mean(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
        itm = itm_result[0];
        itm_unc = itm_result[1];
        h_out->SetBinContent(i, itm);
        h_out->SetBinError(i, itm_unc);
        h_1d_temp->Delete();
    }


    // Clean up
    delete h;
    rfile->Close();

}


void profile_3d_proj(TH2D* h_out, const char* input_file, int dimx, int dimy, int Nx, int Ny, int tpc, int plane, int Ndim) {

    // Open the input file
    TFile* rfile = TFile::Open(input_file);
    if (!rfile || rfile->IsZombie()) {
        std::cerr << "Error opening file: " << input_file << std::endl;
        return;
    }

    // Get the histogram
    int idx = 3 * tpc + plane;
    std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
    const char* cstr = num_str.c_str();
    THnSparseD* h = (THnSparseD*)rfile->Get(cstr);
    if (!h) {
        std::cerr << "Histogram not found: hwidth" << idx << std::endl;
        rfile->Close();
        return;
    }

    // Loop over the bins in the selected dimension
    int dx = h->GetAxis(dimx)->GetNbins() / Nx;
    int dy = h->GetAxis(dimy)->GetNbins() / Ny;

    for (int binx = 1; binx < Nx + 1; ++binx) {
        for (int biny = 1; biny < Ny + 1; ++biny) {

            h->GetAxis(dimx)->SetRange(binx*dx, binx*dx +dx);
            h->GetAxis(dimy)->SetRange(biny*dy, biny*dy + dy);
            TH1D* h_1d_temp = h->Projection(Ndim -1);
            Float_t itm = 0.; // iterative truncated mean (ITM)
            Float_t itm_unc = 0.; // uncertainty of ITM
            Double_t itm_result[2];
            iterative_truncated_mean(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
            itm = itm_result[0];
            itm_unc = itm_result[1];

            for (int i = binx*dx; i < binx*dx +dx + 1; ++i) {
                for (int j = biny*dy; j < biny*dy +dy + 1; ++j) {
                    h_out->SetBinContent(i, j, itm);
                    h_out->SetBinError(i, j, itm_unc);
                }
            }
            
            h_1d_temp->Delete();
            
        }
    }
    
    // Clean up
    delete h;
    rfile->Close();

}

