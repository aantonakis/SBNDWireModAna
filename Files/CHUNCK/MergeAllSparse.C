#include <TFile.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

void MergeAllSparse(const char* filelist = "files.txt",
                    const char* outname = "merged.root")
{
    std::ifstream infile(filelist);
    if (!infile.is_open()) {
        Error("MergeAllSparse", "Cannot open file list: %s", filelist);
        return;
    }

    std::vector<std::string> filenames;
    std::string line;
    while (infile >> line) filenames.push_back(line);
    infile.close();

    if (filenames.empty()) {
        Error("MergeAllSparse", "File list is empty!");
        return;
    }

    // Open the first file to discover all THnSparseD keys
    TFile* f0 = TFile::Open(filenames.front().c_str(), "READ");
    if (!f0 || f0->IsZombie()) {
        Error("MergeAllSparse", "Cannot open first file: %s", filenames.front().c_str());
        return;
    }

    std::map<std::string, THnSparseD*> sums;

    TIter nextkey(f0->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if (cl->InheritsFrom("THnSparseD")) {
            std::string hname = key->GetName();
            THnSparseD* h0 = (THnSparseD*)f0->Get(hname.c_str());
            if (!h0) continue;

            THnSparseD* hsum = (THnSparseD*)h0->Clone(Form("hsum_%s", hname.c_str()));
            //hsum->SetDirectory(nullptr);
            hsum->Reset();
            sums[hname] = hsum;
        }
    }
    f0->Close();

    printf("Discovered %zu THnSparseD histograms in first file.\n", sums.size());

    // Loop over all files and merge
    int nfiles = 0;
    for (const auto& fname : filenames) {
        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            Warning("MergeAllSparse", "Skipping bad file: %s", fname.c_str());
            continue;
        }

        for (auto& [hname, hsum] : sums) {
            THnSparseD* htemp = (THnSparseD*)f->Get(hname.c_str());
            if (!htemp) continue;
            hsum->Add(htemp);
        }

        f->Close();
        ++nfiles;

        if (nfiles % 100 == 0)
            printf("  ... merged %d files\n", nfiles);
    }

    printf("Merged %d files total.\n", nfiles);

    // Write the merged histograms
    TFile* fout = TFile::Open(outname, "RECREATE");
    for (auto& [hname, hsum] : sums)
        if (hsum) hsum->Write(hsum->GetName());
    fout->Close();

    printf("Merged THnSparseD histograms written to %s\n", outname);
}

