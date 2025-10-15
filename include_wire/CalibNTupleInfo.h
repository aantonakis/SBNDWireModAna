#ifndef CALIB_NTUPLE_INFO_H
#define CALIB_NTUPLE_INFO_H

// Standard Library Includes
#include <iostream>
#include <fstream>

// ROOT Includes
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;

struct MyCalib {

    TTreeReader reader; 
    TTreeReaderValue<int> run;
    TTreeReaderValue<int> subrun;
    TTreeReaderValue<int> evt;
    TTreeReaderValue<int> selected;
    TTreeReaderValue<int> whicht0;

    TTreeReaderValue<float> trk_dirx;
    TTreeReaderValue<float> trk_diry;
    TTreeReaderValue<float> trk_dirz;

    TTreeReaderArray<unsigned short> tpc[kNplanes];
    TTreeReaderArray<float> goodness[kNplanes];

    TTreeReaderArray<float> x[kNplanes];
    TTreeReaderArray<float> y[kNplanes];
    TTreeReaderArray<float> z[kNplanes];

    TTreeReaderArray<float> width[kNplanes];
    TTreeReaderArray<float> integral[kNplanes];
    TTreeReaderArray<bool> ontraj[kNplanes];

    TTreeReaderArray<float> dirx[kNplanes];
    TTreeReaderArray<float> diry[kNplanes];
    TTreeReaderArray<float> dirz[kNplanes];

    TTreeReaderArray<float> dqdx[kNplanes];
    TTreeReaderArray<float> rr[kNplanes]; 

    MyCalib(TChain *fChain)
        : reader(fChain),
          run(reader, "meta.run"),
          subrun(reader, "meta.subrun"),
          evt(reader, "meta.evt"), 
	  selected(reader, "trk.selected"),
	  whicht0(reader, "trk.whicht0"),
	  trk_dirx(reader, "trk.dir.x"),
	  trk_diry(reader, "trk.dir.y"),
	  trk_dirz(reader, "trk.dir.z"),


	  // Arrays	
	  tpc{              
	    TTreeReaderArray<unsigned short>(reader, "trk.hits0.h.tpc"),
            TTreeReaderArray<unsigned short>(reader, "trk.hits1.h.tpc"),
            TTreeReaderArray<unsigned short>(reader, "trk.hits2.h.tpc")
	  },
	  
	  goodness{
            TTreeReaderArray<float>(reader, "trk.hits0.h.goodness"),
            TTreeReaderArray<float>(reader, "trk.hits1.h.goodness"),
            TTreeReaderArray<float>(reader, "trk.hits2.h.goodness")
	  },
	  
  	  x{
	    TTreeReaderArray<float>(reader, "trk.hits0.h.sp.x"),
	    TTreeReaderArray<float>(reader, "trk.hits1.h.sp.x"),
	    TTreeReaderArray<float>(reader, "trk.hits2.h.sp.x")
	  },
  	  y{
	    TTreeReaderArray<float>(reader, "trk.hits0.h.sp.y"),
	    TTreeReaderArray<float>(reader, "trk.hits1.h.sp.y"),
	    TTreeReaderArray<float>(reader, "trk.hits2.h.sp.y")
	  },
  	  z{
	    TTreeReaderArray<float>(reader, "trk.hits0.h.sp.z"),
	    TTreeReaderArray<float>(reader, "trk.hits1.h.sp.z"),
	    TTreeReaderArray<float>(reader, "trk.hits2.h.sp.z")
	  },

	  width{
	    TTreeReaderArray<float>(reader, "trk.hits0.h.width"),
	    TTreeReaderArray<float>(reader, "trk.hits1.h.width"),
	    TTreeReaderArray<float>(reader, "trk.hits2.h.width")
	  },

	  integral{
	    TTreeReaderArray<float>(reader, "trk.hits0.h.integral"),
	    TTreeReaderArray<float>(reader, "trk.hits1.h.integral"),
	    TTreeReaderArray<float>(reader, "trk.hits2.h.integral")
	  },
	  
	  ontraj{
	    TTreeReaderArray<bool>(reader, "trk.hits0.ontraj"),
	    TTreeReaderArray<bool>(reader, "trk.hits1.ontraj"),
	    TTreeReaderArray<bool>(reader, "trk.hits2.ontraj"),
	  },
	  
	  dirx{
            TTreeReaderArray<float>(reader, "trk.hits0.dir.x"),
            TTreeReaderArray<float>(reader, "trk.hits1.dir.x"),
            TTreeReaderArray<float>(reader, "trk.hits2.dir.x"),
	  },
	  diry{
            TTreeReaderArray<float>(reader, "trk.hits0.dir.y"),
            TTreeReaderArray<float>(reader, "trk.hits1.dir.y"),
            TTreeReaderArray<float>(reader, "trk.hits2.dir.y"),
	  },
	  dirz{
            TTreeReaderArray<float>(reader, "trk.hits0.dir.z"),
            TTreeReaderArray<float>(reader, "trk.hits1.dir.z"),
            TTreeReaderArray<float>(reader, "trk.hits2.dir.z"),
	  },
	  
	  dqdx{
	    TTreeReaderArray<float>(reader, "trk.hits0.dqdx"),
	    TTreeReaderArray<float>(reader, "trk.hits1.dqdx"),
	    TTreeReaderArray<float>(reader, "trk.hits2.dqdx"),
	  },
	  
	  rr{
	    TTreeReaderArray<float>(reader, "trk.hits0.rr"),
	    TTreeReaderArray<float>(reader, "trk.hits1.rr"),
	    TTreeReaderArray<float>(reader, "trk.hits2.rr")
	  }

	  {}
};

#endif

