
#ifndef PROFILE_H
#define PROFILE_H

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
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
    int profileDim,
    int projDim,
    int minCounts = 1000)
{
    ProfileResult result;
    if (!h) {
        std::cerr << "Null THnSparseD provided!\n";
        return result;
    }

    TAxis* axis = h->GetAxis(profileDim);
    int nBins = axis->GetNbins();
    std::vector<double> binEdges;   // store dynamic edges
    std::vector<double> binCounts;  // store counts per dynamic bin

    binEdges.push_back(axis->GetBinLowEdge(1)); // first edge

    int startBin = 1;
    while (startBin <= nBins) {
        int endBin = startBin;
        double totalEntries = 0.0;

        // expand range until minCounts satisfied or we run out of bins
        
        while (endBin <= nBins && totalEntries < minCounts) {
            totalEntries += axis->GetBinContent(endBin);
            endBin++;
        }

	// set range for this slice
	axis->SetRange(startBin, endBin - 1);

	// Project
        TH1D* proj = h->Projection(projDim);
        proj->SetDirectory(nullptr);
        proj->SetName(Form("proj_dim%d_%d_%d", projDim, startBin, endBin - 1));
        result.projections.push_back(proj);

	// Record binning info
        binEdges.push_back(axis->GetBinUpEdge(endBin - 1));
        binCounts.push_back(totalEntries);

	startBin = endBin;

    }
    // reset axis
    axis->SetRange();

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



void SaveProfileResult1D(const ProfileResult& result, const std::string& outFile) {
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

void SaveAllProfileResult1D(const std::map<std::string, ProfileResult>& results,
                           const std::string& outFile) {
    TFile f(outFile.c_str(), "RECREATE");

    for (const auto& [name, res] : results) {
        TDirectory* resDir = f.mkdir(name.c_str());
        SaveProfileResultToDir(res, resDir);
    }

    f.Close();
}

#endif
