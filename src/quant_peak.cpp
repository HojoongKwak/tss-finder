#include "bg.h"					// Custom reading bedgraph files
#include "datafile.h"			// Custom reading tab delimited lists
#include "arguments.h"			// Parsing command line arguments
#include "smooth.h"				// Custom Gaussian smoothing functions
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <thread>
#include <functional>
using namespace std;

// Declarations
struct args		// Arguments
{
	char *pfn,*mfn;			// Bedgraph filenames
	char *gfn,*ofn;			// Output filename
	int nth;
	int res,sh,bin;				// resolution
	// Initialization
	args() { nth=6; res=300; sh=80; bin=10; };
	bool get(int argc, char *argv[]);	// Function to parse arguments
};

// Gene coordinates format
struct gene_coordinates
{
	vector<string> chrom;		// Chromosome
	vector<string> strand;		// Strand [+/-]
	vector<int> start;			// Start positions
	vector<int> end;			// End positions
	int n;						// Number of genes
	void load(char *fn);		// Read from a tab delimited text
};

// Functions used in main
void get_data(vector<vector<float> > &pden,vector<vector<float> > &mden,bgdata &p,bgdata &m,args &a,int chrgroup);
void get_peak(vector<vector<float> > &pden,vector<vector<float> > &mden, vector<vector<int> > &peakpos,
	vector<vector<float> > &psum, vector<vector<float> > &msum, args &a, int chrgroup);

int main(int argc,char *argv[])
{
	// Parsing arguments
	args a;
	if(!a.get(argc,argv)) return 0;

	gene_coordinates g;
	g.load(a.gfn);	
	// Load bedgraph data
	bgdata p,m;
	p.load(a.pfn);		// Plus strand data
	m.load(a.mfn);	// Minus strand data if provided
	int nchr=p.chr.size();
		
	// Arrays of bedgraph data values for every gene
	vector<vector<float> > pden(nchr),mden(nchr);	// Data arrays for bedgraph values
	// Read data into data array, retrieve read counts
	// Use multithreading
	vector<thread> thr(a.nth);
	cout<<"Reading bedgraphs"<<endl;
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_data,ref(pden),ref(mden),ref(p),ref(m),ref(a),i);
	for(int i=0;i<a.nth;i++)
	{
		thr[i].join();
	}
	cout<<"Finished reading bedgraphs"<<endl;

	vector<vector<int> > peakpos(nchr);
	for(int i=0;i<g.n;++i)
	{
		int c=p.getChrID(g.chrom[i]);
		if(c>=0&&c<nchr)
		{
			if(g.strand[i]=="+") peakpos[c].push_back(g.start[i]/a.bin);
			else peakpos[c].push_back(g.end[i]/a.bin);
		}
	}

	vector<vector<float> > psum(nchr);
	vector<vector<float> > msum(nchr);
	cout<<"Scanning peaks"<<endl;
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_peak,ref(pden),ref(mden),ref(peakpos),ref(psum),ref(msum),ref(a),i);
	for(int i=0;i<a.nth;i++)
	{
		thr[i].join();
	}
	cout<<"Finished scanning peaks"<<endl;

	ofstream out(a.ofn);
	out<<"chr\tpos\tplus\tminus"<<endl;
	for(int i=0;i<peakpos.size();i++)
	{
		for(int j=0;j<peakpos[i].size();j++)
			out<<p.chr[i]<<"\t"<<peakpos[i][j]*a.bin<<"\t"<<psum[i][j]<<"\t"<<msum[i][j]<<endl;
	}
}
	   
// Parsing arguments
bool args::get(int argc, char *argv[])
{
	// Parsing arguments
	arg::push("-p","plus bedgraph filename",pfn);
	arg::push("-m","minus bedgraph filename (optional)",mfn);
	arg::push("-g","reference list filename",gfn);
	arg::push("-o","output filename",ofn);
	arg::push("-r","peak resolution (default=300)",res,true);
	arg::push("-b","bin size (default=10)",bin,true);
	arg::push("-sh","shift (default=80)",sh,true);
	arg::push("-nth","Number of threads (default=6)",nth,true);
	return(arg::get(argc,argv));
};

// Functions for loading gene_coordinates
void gene_coordinates::load(char *fn)
{
	datafile d;
	d.load(fn);
	d.get_data(chrom,0);
	d.get_data(strand,1);
	d.get_data(start,2);
	d.get_data(end,3);
	n=chrom.size();
};

// Function for making a data array from begraph files
void get_data(vector<vector<float> > &pden,vector<vector<float> > &mden,bgdata &p,bgdata &m,args &a,int chrgroup)
{
	int sw=a.res/a.bin;
	smoothFunction sf(sw);
	for(int chr=0;chr<p.chr.size();chr++)		// Go through all chromosomes
	{
		if((chr%a.nth)!=chrgroup) continue;	// For multithreading
		string chrname=p.chr[chr];
		p.getlevel(pden[chr],chrname,a.sh/a.bin,p.end[p.getChrID(chrname)].back(),a.bin);
		if(m.getChrID(chrname)!=-1)
			m.getlevel(mden[chr],chrname,-a.sh/a.bin,m.end[m.getChrID(chrname)].back(),a.bin);
		int chrlen=pden[chr].size();
		if(chrlen<mden[chr].size()) chrlen=mden[chr].size();
		if(chrlen<sw) chrlen=sw;
		pden[chr].resize(chrlen);
		mden[chr].resize(chrlen);
		sf.apply(pden[chr]);
		sf.apply(mden[chr]);
	}
};


void get_peak(vector<vector<float> > &pden,vector<vector<float> > &mden, vector<vector<int> > &peakpos,
	vector<vector<float> > &psum, vector<vector<float> > &msum, args &a,int chrgroup)
{
	int nchr=pden.size();
	for(int chr=0;chr<nchr;++chr)
	{
		if((chr%a.nth)!=chrgroup) continue;
		for(int i=0;i<peakpos[chr].size();i++)
		{
			int pos=peakpos[chr][i];
			if(pos>=0&&pos<pden[chr].size())
			{
				psum[chr].push_back(pden[chr][pos]);
				msum[chr].push_back(mden[chr][pos]);
			}
			else
			{
				psum[chr].push_back(0);
				msum[chr].push_back(0);
			}
		}
	}
};
