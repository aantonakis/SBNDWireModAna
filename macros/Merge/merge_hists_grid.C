#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "Math/Vector3D.h"

#include "CalibrationStandard.h"
#include "CalibNTupleInfo.h"
#include "Angles.h"

#include "SelectionWire.h"

using ROOT::Math::XYZVector;

std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");

void printMemoryUsage()
{
    ProcInfo_t info;
    gSystem->GetProcInfo(&info);

    std::cout << "--------------------------------------------------\n";
    std::cout << "Memory usage (current process):\n";
    std::cout << "  Resident set size (RSS): " << info.fMemResident << " MB\n";
    std::cout << "  Virtual memory size   : " << info.fMemVirtual  << " MB\n";
    std::cout << "  CPU user time         : " << info.fCpuUser << " s\n";
    std::cout << "  CPU sys time          : " << info.fCpuSys  << " s\n";
    std::cout << "--------------------------------------------------\n";
}


void merge_hists_grid(TString list_file, TString out_suffix,

  std::vector<int> dim = {0} // dimesnions to project 

) {

    // File List Management
    TString input_file_dir = getenv("DATA_PATH");
    TString sample_list_dir = getenv("SAMPLE_PATH");
    TString sample_list_label = getenv("FILELIST_LABEL");
    
    TString fileListPath = sample_list_dir + "/" + list_file;
    cout << "Opening : " << fileListPath << endl;

    std::ifstream file(fileListPath.Data());  // Convert TString to const char*
    if (!file) {
      cout << "File does not exist: " << fileListPath << endl;
      cout << "Exiting" << endl;
      return;
    }

    std::vector<string> files;
    
    ifstream in(fileListPath);
    string fileName;
    while(getline(in, fileName)) {
      //std::cout << "Adding file : " << fileName.c_str() << std::endl;
      files.push_back(fileName);
    }
    
    // initialize histograms
    THnSparseD* h[kNplanes * kNTPCs];
    THnSparseD* hTracks[kNplanes * kNTPCs];

    std::cout << "Finding first valid histograms ..." << std::endl;
    int start = 0;
    bool stop = false;
    for (int i = 0; i < files.size(); ++i) {
      TFile* f = TFile::Open(files[i].c_str(), "READ");
      if (!f) {
        start += 1;
        continue;
      } 
      THnSparseD* h_temp = (THnSparseD*)f->Get("hHit0");
      if (!h_temp) {
        h_temp->Delete();
        start += 1;
        continue;
      }  
      else {
        h_temp->Delete();
        std::cout << "Found valid first hist --> Adding all starting hists" << std::endl;
        stop = true;
        for (unsigned j = 0; j < kNplanes * kNTPCs; j++) {
          THnSparseD* h_temp = (THnSparseD*)f->Get(Form("hHit%d", j));
          THnSparseD* h_temp_trk = (THnSparseD*)f->Get(Form("hTrack%d", j));
          h[j] = static_cast<THnSparseD*>( h_temp->Projection(dim.size(), dim.data()) );
          hTracks[j] = static_cast<THnSparseD*>( h_temp_trk->Projection(dim.size(), dim.data()) );
          h[j]->SetName(Form("hHit%d", j));
          hTracks[j]->SetName(Form("hTrack%d", j));
          h_temp->Delete();
          h_temp_trk->Delete();
        } 
      }
      f->Close();
      delete f;
      if (stop) break;
    }

    std::cout << "Adding remianing hists from the " << files.size() << " files ..." << std::endl;
    // loop over the remainder of the files and add
    for (int i = start; i < files.size(); ++i) {
      TFile* f = TFile::Open(files[i].c_str(), "READ");
      if (!f) {
        continue;
      } 
      std::cout << "Adding file " << i << std::endl;
      if (i % 10 == 0) printMemoryUsage();

      for (unsigned j = 0; j < kNplanes * kNTPCs; j++) {
        THnSparseD* h_temp = (THnSparseD*)f->Get(Form("hHit%d", j));
        THnSparseD* h_temp_trk = (THnSparseD*)f->Get(Form("hTrack%d", j));
        if ((!h_temp) || (!h_temp_trk)) continue;
        h[j]->Add(static_cast<THnSparseD*>( h_temp->Projection(dim.size(), dim.data()) ));
        hTracks[j]->Add(static_cast<THnSparseD*>( h_temp_trk->Projection(dim.size(), dim.data()) ));
        h_temp->Delete();
        h_temp_trk->Delete();
      }
      f->Close();
      delete f;
    }


    std::cout << "About to write histograms to the output file" << std::endl;

    TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
    TString output_file_name = output_rootfile_dir + "/output_merge_hists_" + out_suffix + ".root";
    out_rootfile = new TFile(output_file_name, "RECREATE");
    out_rootfile -> cd();
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
	std::cout << "Writing histograms for plane " << i << std::endl;
        h[i]->Write();
        hTracks[i]->Write();
    }
   
    out_rootfile->Close();

} // end of main


std::vector<TString> filenames_from_input(const TString& input_arg, int nmax=-1) {
    std::vector<TString> filenames;
    if (input_arg.EndsWith(".root")) {
        filenames.push_back(input_arg);
        return filenames;
    }

    // read from file
    size_t nfiles = 0;
    std::ifstream ifile(input_arg);
    std::string line;
    while (std::getline(ifile, line)) {
        nfiles++;
        fprintf(stdout, "Adding file %zu: %s...\n", nfiles, line.c_str());
        filenames.push_back(TString(line));
        if (nfiles >= nmax && nmax > 0) break;
    }
    return filenames;
}


TString basename_prefix(const TString& input, const TString& prefix, const TString& suffix) {
    // remove path from filename and return new string with prefix or suffix added before extension
    TPRegexp re(".*/(.*)");
    TObjArray* matches = re.MatchS(input);
    TString result((static_cast<TObjString*>(matches->At(1)))->String());
    matches->Delete();
    return prefix + result;
}




