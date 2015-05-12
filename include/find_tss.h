#include "bed.h"
#include "bg.h"

#include <string>
#include <vector>

using namespace std;

typedef vector<vector<float> > float_2d;
typedef vector<vector<int> > int_2d;

// Parameter definition
struct par
{
	char *pfn,*mfn;			// Bedgraph filenames
	string ofn; 			// Output file prefix
	int nth;                // Threads to use (6)
	int window;             // TSS window size (20 bp)
    int shStart,shEnd;      // Range of shift distances to consider (-120 to +120)
    int bin;			    // Bin sizes (5)
    int cutoff;             // Cutoff readcount per bin

    par();                  // Constructor with parameter initializations
	bool get(int argc, char *argv[]);	// Function to parse arguments
};

// Functions used in main
void get_bgdata(float_2d &pden, float_2d &mden, bgdata &p, bgdata &m, par &a);
void make_shbed(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, vector<string> &chromName, par &a);
void save_bed(vector<vector<bedTrack> > &den, string outputPrefix); 
