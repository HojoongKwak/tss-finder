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
  string pfn, mfn;
  string ofn; 			// Output file prefix
  int nth;                // Threads to use (6)
  int window;             // TSS window size (20 bp)
  int shStart,shEnd;      // Range of shift distances to consider (-120 to +120)
  int bin;			    // Bin sizes (5)

  par();                  // Constructor with parameter initializations
  bool check();
};

// Functions used in main
void get_bgdata(float_2d &pden, float_2d &mden, bgdata &p, bgdata &m, par &a);
void make_shbed(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, par &a, int shift);

