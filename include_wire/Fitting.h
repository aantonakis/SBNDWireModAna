// Header File for itm_fit.cc
#ifndef FITTING_H
#define FITTING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>

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

// Langau Fit Reference: https://root.cern/doc/master/langaus_8C.html

// ---------------- LANGAU FUNCTIONS ---------------------- //

Double_t langaufun(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.398942280401;// Control constants
    // Double_t mpshift = -0.22278298l
    Double_t np = 500.0;
    Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    mpc=par[1];
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow)/np;

    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}


TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
 
   int i;

   //char FunName[100];
   //sprintf(FunName,"Fitfcn_%s",his->GetName());

   std::string FunNameS = std::string("Fitfcn_") + his->GetName();
    const char* FunName = FunNameS.c_str();
 
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;
 
   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
 
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }
 
   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
 
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf
 
   return (ffit);              // return fit function
 
}

int langaupro(double *params, double &maxx, double &FWHM) {
 
   // Searches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.
 
   double p,x,fy,fxr,fxl;
   double step;
   double l,lold;
   int i = 0;
   int MAXCALLS = 10000;
 
 
   // Search for maximum
 
   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-1);
 
   maxx = x;
 
   fy = l/2;
 
 
   // Search for right x location of fy
 
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-2);
 
   fxr = x;
 
 
   // Search for left x location of fy
 
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-3);
 
 
   fxl = x;
 
   FWHM = fxr - fxl;
   return (0);
}


// ------------------------------------------------------ //


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


void hist_mean_std_error(const TH1* h, Double_t* result) {
    //Int_t nbins = h->GetNbinsX();
    if (h->Integral() == 0) return; 

    Double_t mean = h->GetMean();
    Double_t sigma = h->GetRMS();
    Int_t N = h->GetEntries();
    Int_t NI = h->Integral();
    //if (NI != N) {
    //    std::cout << "Entries and Integral Don't Match?" << std::endl;
    //    std::cout << "Entries: " << N << " Integral: " << NI << std::endl;
    //}
    Double_t err = sigma / sqrt(NI);

    result[0] = mean;
    result[1] = err;
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


void iterative_truncated_mean_std_err(TH1* h, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result) {
    if (sig_down > sig_up) {
        fprintf(stderr, "Warning: reversing iterative truncated mean limits [%.2e,%.2e]\n", sig_down, sig_up);
        Double_t tmp = sig_down;
        sig_down = sig_up;
        sig_up = tmp;
    }

    Double_t mean = h->GetMean();
    Double_t sd = h->GetRMS();
    //std::cout << mean << ", " << sd << ", " << result[0] << "\n";

    // return mean, std err of iterative truncated mean
    if (TMath::Abs(result[0] - mean) < tol) {
        hist_mean_std_error(h, result);
        return;
    }
    result[0] = mean;
    result[1] = sd;
    
    // get positions to cut based on quantiles
    Double_t probs[1] = { 0.5 };
    Double_t median[1];
    h->GetQuantiles(1, median, probs);

    //std::cout << "median: " << median[0] << "\n";
    //std::cout << median[0] + sig_down * sd << ", " << median[0] + sig_up * sd << "\n";
    
    TH1* hnew = (TH1*)h->Clone();
    hnew->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinLowEdge(i) > median[0] + sig_up * sd || h->GetBinLowEdge(i) + h->GetBinWidth(i) < median[0] + sig_down * sd) continue;
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }
    
    iterative_truncated_mean_std_err(hnew, sig_down, sig_up, tol, result);
    delete hnew;
    return;
}


void FitPoly3(TH1* h, Double_t* result) {
	if (result[0] == 0) return;
	
	int firstBin = h->FindFirstBinAbove(0);
	int lastBin  = h->FindLastBinAbove(0);
	double xmin = h->GetBinLowEdge(firstBin);
	double xmax = h->GetBinLowEdge(lastBin) + h->GetBinWidth(lastBin);
	
	TF1* f = new TF1("f", "pol3", xmin, xmax);
    	TFitResultPtr r = h->Fit(f, "SR");

	if (r->Status() != 0) {
          std::cout << "Fit failed with status = " << r->Status() << "\n";
	  return;
	} else {
          std::cout << "Fit succeeded!\n";
        }
	float goodness = r->Chi2() / r->Ndf();
	if (goodness > 3) {
	  std::cout << "BAD FIT chi/N = " << goodness << std::endl;
	  return;
	}

        double p1 = f->GetParameter(1);
        double p2 = f->GetParameter(2);
        double p3 = f->GetParameter(3);

        double disc = 4*p2*p2 - 12*p3*p1;  // discriminant
        if (disc < 0) {
          std::cout << "No real roots of derivative!\n";
          return;
        }

        double root1 = (-2*p2 + std::sqrt(disc)) / (6*p3);
        double root2 = (-2*p2 - std::sqrt(disc)) / (6*p3);

        auto secondDeriv = [&](double x) {
          return 2*p2 + 6*p3*x;
        };

        double peakX = NAN;
        if (secondDeriv(root1) < 0) { 
		peakX = root1;
	}
        else if (secondDeriv(root2) < 0) {
		 peakX = root2;
	}
	else {
		std::cout << "BAD 2nd Deriv" << std::endl;
		return;
	}
        std::cout << "Peak at x = " << peakX << std::endl;
	std::cout << "Peak from ITM " << result[0] << std::endl;
    	//double root = TMath::BrentRoot([&](double x){ return deriv(x); }, xmin, xmax);
	if (peakX == NAN) {
	  std::cout << "PEAK IS NAN" << std::endl;
	  return;
	}
	
	result[0] = peakX;
}

void iterative_truncated_mean_poly3(TH1* h, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result) {
    if (sig_down > sig_up) {
        fprintf(stderr, "Warning: reversing iterative truncated mean limits [%.2e,%.2e]\n", sig_down, sig_up);
        Double_t tmp = sig_down;
        sig_down = sig_up;
        sig_up = tmp;
    }

    Double_t mean = h->GetMean();
    Double_t sd = h->GetRMS();
    //std::cout << mean << ", " << sd << ", " << result[0] << "\n";

    // return mean, std err of iterative truncated mean
    if (TMath::Abs(result[0] - mean) < tol) {
        hist_mean_std_error(h, result);
	if (result[0] < tol) return;
	FitPoly3(h, result);
        return;
    }
    result[0] = mean;
    result[1] = sd;
    
    // get positions to cut based on quantiles
    Double_t probs[1] = { 0.5 };
    Double_t median[1];
    h->GetQuantiles(1, median, probs);

    //std::cout << "median: " << median[0] << "\n";
    //std::cout << median[0] + sig_down * sd << ", " << median[0] + sig_up * sd << "\n";
    
    TH1* hnew = (TH1*)h->Clone();
    hnew->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinLowEdge(i) > median[0] + sig_up * sd || h->GetBinLowEdge(i) + h->GetBinWidth(i) < median[0] + sig_down * sd) continue;
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }
    
    iterative_truncated_mean_poly3(hnew, sig_down, sig_up, tol, result);
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

void profile_2d_proj_std_err(TH1D* h_out, const char* input_file, int dim, int tpc, int plane, int Ndim) {
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
        iterative_truncated_mean_std_err(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
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

void profile_2d_proj_poly3(TH1D* h_out, const char* input_file, int dim, int tpc, int plane, int Ndim) {
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
        //iterative_truncated_mean_std_err(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
        iterative_truncated_mean_poly3(h_1d_temp, -2, 1.75, 1.0e-4, itm_result);
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

#endif
