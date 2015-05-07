#include "bg.h"					// Custom reading bedgraph files
#include "params.h"
#include "smooth.h"				// Custom Gaussian smoothing functions
#include "bed.h"
#include "Error.h"

#include "find_tss_hk.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <thread>
#include <functional>

using namespace std;

int main(int argc,char *argv[])
{
  par pa;  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required Options",NULL)
    LONG_STRING_PARAM("pl",&pa.pfn,"Plus strand bedgraph filename")
    LONG_STRING_PARAM("mn",&pa.mfn,"Minus strand bedgraph filename")
    LONG_STRING_PARAM("out",&pa.ofn,"Output file prefix")
    
    LONG_PARAM_GROUP("Additional Options",NULL)
    LONG_INT_PARAM("win",&pa.window,"TSS window size")
    LONG_INT_PARAM("min-shift",&pa.shStart,"Maximum convergent (or minimum divergent) shift range")
    LONG_INT_PARAM("max-shift",&pa.shEnd,"Maximum divergent shift range")
    LONG_INT_PARAM("bin",&pa.bin,"Bin size")
    LONG_INT_PARAM("threads",&pa.nth,"Number of threads")        
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( !pa.check() ) {
    error("Missing required parameter(s)");
  }  
  
  // Load bedgraph data
  bgdata bp,bm;
  bp.load((char*)pa.pfn.c_str());                  // Plus strand data
  if(!pa.mfn.empty()) bm.load((char*)pa.mfn.c_str());	    // Minus strand data if provided
  else bm=bp;
  int nchr=bp.chr.size();
  
  // Arrays of bedgraph data values for every chromosome
  float_2d pden(nchr),mden(nchr);	// Data arrays for bedgraph values
  
  // Read read counts into data array
  cerr<<"Reading bedgraph data..."<<endl;
  get_bgdata(pden,mden,bp,bm,pa);     // Sum of read counts within the window size from binned positions
  
  // Repeat routine for the whole shift range
  for(int i=pa.shStart;i<=pa.shEnd;i+=pa.bin)
    {
        cerr<<"Analyzing shift size = "<<i<<" bp..."<<endl;
        vector<vector<bedTrack> > den(nchr);
        make_shbed(den, pden, mden, pa, i);
    }
};


    
