#include "bg.h"					// Custom reading bedgraph files
#include "arguments.h"			// Parsing command line arguments
#include "smooth.h"				// Custom Gaussian smoothing functions
#include "bed.h"

#include "find_tss.h"

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
	// Parsing arguments to parameters
	par pa;
	if(!pa.get(argc,argv)) return 0;
	
	// Load bedgraph data
	bgdata bp,bm;
	bp.load(pa.pfn);                    // Plus strand data
	if(pa.mfn) bm.load(pa.mfn);	        // Minus strand data if provided
	else bm=bp;
    int nchr=bp.chr.size();
    vector<string> chromName=bp.chr;    // Chromosome names
		
	// Arrays of bedgraph data values for every chromosome
	float_2d pden(nchr),mden(nchr);	    // Data arrays for bedgraph values
	
    // Read read counts into data array
	cerr<<"Reading bedgraph data..."<<endl;
    get_bgdata(pden,mden,bp,bm,pa);     // Sum of read counts within the window size from binned positions
    
	// Analyze TSS peaks
    cerr<<"Analyzing TSS peak pairs..."<<endl;
    vector<vector<bedTrack> > den(nchr);
    make_shbed(den, pden, mden, chromName, pa);     // Make bed tracks
    save_bed(den, pa.ofn);                          // Save bed tracks
};


    
