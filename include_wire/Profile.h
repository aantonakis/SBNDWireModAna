
#ifndef PROFILE_H
#define PROFILE_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVector;


struct ProfileResult1D {
    std::vector<TH1D*> projections;  // the 1D projections over projDim
    TH1D* summary;                   // dynamic binning summary histogram
};

struct ProfileResult2D {
    std::vector<TH2D*> projections;  // the 2D projections over projDim
    TH2D* summary;                   // dynamic binning summary histogram
};


/// Profile a THnSparseD over one dimension with dynamic binning.
///// Also returns a summary TH1D with non-uniform binning and counts per bin.
///// @param h          Input THnSparseD
///// @param profileDim Dimension index to profile over
///// @param projDim    Dimension index to project into TH1D
///// @param minCounts  Minimum counts required per bin (default 1000)
///// @return ProfileResult with projections and summary histogram

ProfileResult1D ProfileSparseDynamic1D(
    THnSparseD* h,
    THnSparseD* hT, // track count histogram
    int profileDim,
    int projDim,
    int minCounts = 1000)
{
    ProfileResult1D result;
    if (!h) {
        std::cerr << "Null THnSparseD provided!\n";
        return result;
    }
    TH1D* h_1d_axis = hT->Projection(profileDim);
    std::cout << "Got the 1D axis projection" << std::endl;

    //TAxis* axis = h->GetAxis(profileDim);
    int nBins = h_1d_axis->GetNbinsX();
    std::vector<double> binEdges;   // store dynamic edges
    std::vector<double> binCounts;  // store counts per dynamic bin

    binEdges.push_back(h_1d_axis->GetXaxis()->GetBinLowEdge(1)); // first edge
    std::cout << "Start dynamic binning with " << nBins << " bins, min counts: " << minCounts << std::endl;
    int startBin = 1;
    int newBin = 1;
    while (startBin <= nBins) {
        int endBin = startBin;
        double totalEntries = 0.0;
        std::cout << "Starting new dynamic bin at bin " << startBin << std::endl;
        // expand range until minCounts satisfied or we run out of bins
        while (endBin <= nBins && totalEntries < minCounts) {
            if (h_1d_axis->GetBinContent(endBin) == 0) {
                endBin++;
                break;
            }
            totalEntries += h_1d_axis->GetBinContent(endBin);
            endBin++;
        }
        std::cout << "  Created bin from " << startBin << " to " << endBin-1
                  << " with total counts: " << totalEntries << std::endl;

	    // set range for this slice
	    //axis->SetRange(startBin, endBin - 1);

        // set range to current dynamical bin
        h->GetAxis(profileDim)->SetRange(startBin, endBin-1);
        std::cout << "  Set range on profile dim " << profileDim << std::endl;
	    // Project
        TH1D* proj = h->Projection(projDim);
        std::cout << "  Made projection on dim " << projDim << std::endl;
        proj->SetDirectory(nullptr);
        proj->SetName(Form("proj_%d", newBin));
        result.projections.push_back(proj);
        std::cout << "  Stored projection" << std::endl;
	    // Record binning info
        binEdges.push_back(h_1d_axis->GetXaxis()->GetBinUpEdge(endBin-1));
        binCounts.push_back(totalEntries);

	    startBin = endBin;
        std::cout << "  Moving to next bin at " << startBin << std::endl;
        std::cout << "------------------------" << std::endl;
        newBin += 1;
    }
    // reset axis
    //axis->SetRange();
    h->GetAxis(profileDim)->SetRange();

    // make summary histogram with variable binning
    result.summary = new TH1D("profile_summary",
                              Form("Profile summary over dim %d", profileDim),
                              binEdges.size() - 1,
                              binEdges.data());
    result.summary->SetDirectory(nullptr);
    for (size_t i = 0; i < binCounts.size(); i++) {
        result.summary->SetBinContent(i + 1, binCounts[i]);
    }

    return result;
}


void SaveProfileResult1D(const ProfileResult1D& result, const std::string& outFile) {
    TFile f(outFile.c_str(), "RECREATE");

    result.summary->Write("summary");

    TDirectory* projDir = f.mkdir("projections");
    projDir->cd();
    for (size_t i = 0; i < result.projections.size(); i++) {
        result.projections[i]->Write(); // uses the name we gave earlier
    }

    f.Close();
}

// Save a single ProfileResult into an existing directory
void SaveProfileResult1DToDir(const ProfileResult1D& result, TDirectory* dir) {
  dir->cd();

  // Write summary
  result.summary->Write("summary");

  // Create subdir for projections
  TDirectory* projDir = dir->mkdir("projections");
  projDir->cd();
  for (size_t i = 0; i < result.projections.size(); i++) {
    result.projections[i]->Write();  // already has name like proj_dimX_binY
  }
}

void SaveAllProfileResult1D(const std::map<std::string, ProfileResult1D>& results,
                           const char* outFile) {
    TFile f(outFile, "RECREATE");

    for (const auto& [name, res] : results) {
        TDirectory* resDir = f.mkdir(name.c_str());
        SaveProfileResult1DToDir(res, resDir);
    }

    f.Close();
}

#endif
