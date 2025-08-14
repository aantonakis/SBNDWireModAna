

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


void tpc_track_scraper(const char* input_file, const char* output_file) {

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TNtuple* track_tree = new TNtuple("track_tree", "track_tree",
        "run:subrun:evt:idx:sel:sx:sy:sz:ex:ey:ez:dirx:diry:dirz:rr");
        

    TTreeReader reader("caloskim/TrackCaloSkim", f);
    TTreeReaderValue<int> selected(reader, "trk.selected");
    TTreeReaderValue<int> whicht0(reader, "trk.whicht0");
        
    TTreeReaderValue<int> run(reader, "meta.run");
    TTreeReaderValue<int> evt(reader, "meta.evt");
    TTreeReaderValue<int> subrun(reader, "meta.subrun");

    // track length cut: use residual range on collection only
    TTreeReaderArray<float> rr2(reader, "trk.hits2.rr");

    TTreeReaderValue<float> trk_dirx(reader, "trk.dir.x");
    TTreeReaderValue<float> trk_diry(reader, "trk.dir.y");
    TTreeReaderValue<float> trk_dirz(reader, "trk.dir.z");

    TTreeReaderValue<float> trk_sx(reader, "trk.start.x");
    TTreeReaderValue<float> trk_sy(reader, "trk.start.y");
    TTreeReaderValue<float> trk_sz(reader, "trk.start.z");

    TTreeReaderValue<float> trk_ex(reader, "trk.end.x");
    TTreeReaderValue<float> trk_ey(reader, "trk.end.y");
    TTreeReaderValue<float> trk_ez(reader, "trk.end.z");

    int track_idx = 0;
    while (reader.Next()) {
           
        //if (*selected < 1) continue;

        // For CRT T0 study or TPC T0 study
        if (*whicht0 != 0) continue;

        // skip short tracks
        size_t nhits = rr2.GetSize();
        if (nhits == 0) {
            fprintf(stderr, "Warning: Selected track (idx=%d, selected=%d) with no hits? Run=%d, Subrun=%d, Evt=%d. Skipping!\n", track_idx, *selected, *run, *subrun, *evt);
            continue;
        }
        //if (rr2[nhits - 1] < kTrackCut) continue;

        track_tree->Fill(*run, *subrun, *evt, track_idx, *selected,
            *trk_sx, *trk_sy, *trk_sz,
            *trk_ex, *trk_ey, *trk_ez,
            *trk_dirx, *trk_diry, *trk_dirz,
            rr2[nhits - 1]);
    
         
        track_idx++;
    }

    
    outfile->cd();   
    track_tree->SetDirectory(0);
    track_tree->Write();
    outfile->Close();
    delete track_tree;
    delete f;
    printf("Finished processing tracks.\n");
}


