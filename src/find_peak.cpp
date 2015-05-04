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
	char *gfn,*ofn;			// Gene list and output filenamse
	int mar5,mar3;			// Margins
	int nth;
	int sm;
	bool div;
	// Initialization
	args() { mfn=NULL; mar5=1000; mar3=1000; nth=6; sm=0; div=false; };
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
void get_data(vector<vector<float> > &pden, vector<vector<float> > &mden, gene_coordinates &g
	,bgdata &p,bgdata &m,args &a,int chrgroup);
void smooth_density(vector<vector<float> > &den, smoothFunction &sf,int nth,int genegroup);

int main(int argc,char *argv[])
{
	// Parsing arguments
	args a;
	if(!a.get(argc,argv)) return 0;
	
	// Load bedgraph data
	bgdata p,m;
	p.load(a.pfn);		// Plus strand data
	if(a.mfn) m.load(a.mfn);	// Minus strand data if provided
	
	// Read gene coordinates list	
	gene_coordinates g;
   	g.load(a.gfn);			// Load gene list
	
	smoothFunction sf(a.sm);
	// Arrays of bedgraph data values for every gene
	vector<vector<float> > pden(g.n),mden(g.n);	// Data arrays for bedgraph values
	for(int i=0;i<g.n;i++)
	{
		pden[i].resize(a.mar3+a.mar5);
		mden[i].resize(a.mar3+a.mar5);
	}	
	// Read data into data array, retrieve read counts
	// Use multithreading
	vector<thread> thr(a.nth);
	cout<<"Reading bedgraphs"<<endl;
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_data,ref(pden),ref(mden),ref(g),ref(p),ref(m),ref(a),i);
	for(int i=0;i<a.nth;i++)
	{
		thr[i].join();
	}
	cout<<"Finished reading bedgraphs"<<endl;

	// Smooth density
	if(a.div) 
		for(int i=0;i<a.nth;i++) thr[i]=thread(smooth_density,ref(mden),ref(sf),a.nth,i);
	else
		for(int i=0;i<a.nth;i++) thr[i]=thread(smooth_density,ref(pden),ref(sf),a.nth,i);
	for(int i=0;i<a.nth;i++) thr[i].join();
			
	// output
	ofstream out(a.ofn);
	out<<"chr\tstrand\tstart\tend"<<endl;
	for(int i=0;i<g.n;++i)
	if(!a.div)
	{
		int maxpos=max_element(pden[i].begin(),pden[i].end())-pden[i].begin()-a.mar5;
		if(pden[i][maxpos+a.mar5]==0) maxpos=0;
		if(g.strand[i]=="+") out<<g.chrom[i]<<"\t+\t"<<g.start[i]+maxpos<<"\t"<<g.end[i]<<endl;
		else out<<g.chrom[i]<<"\t-\t"<<g.start[i]<<"\t"<<g.end[i]-maxpos<<endl;
	}
	else
	{
		int minpos=min_element(mden[i].begin(),mden[i].end())-mden[i].begin()-a.mar5;
		if(mden[i][minpos+a.mar5]==0) minpos=0;
		if(g.strand[i]=="+") out<<g.chrom[i]<<"\t+\t"<<g.start[i]+minpos<<"\t"<<g.end[i]<<endl;
		else out<<g.chrom[i]<<"\t-\t"<<g.start[i]<<"\t"<<g.end[i]-minpos<<endl;
	}
}
	   
// Parsing arguments
bool args::get(int argc, char *argv[])
{
	// Parsing arguments
	arg::push("-p","plus bedgraph filename",pfn);
	arg::push("-m","minus bedgraph filename (optional)",mfn,true);
	arg::push("-g","gene list filename",gfn);
	arg::push("-o","output filename",ofn);
	arg::push("-m5","5' margin (default=1000)",mar5,true);
	arg::push("-m3","3' margin (default=1000)",mar3,true);
	arg::push("-sm","smoothing size (default=0)",sm,true);
	arg::push("-div","find divergent (default=false)",div,true);
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
void get_data(vector<vector<float> > &pden,vector<vector<float> > &mden, gene_coordinates &g
	,bgdata &p,bgdata &m,args &a,int chrgroup)
{
	vector<float> pl,mn;
	for(int chr=0;chr<p.chr.size();chr++)		// Go through all chromosomes
	{
		if((chr%a.nth)!=chrgroup) continue;	// For multithreading
		string chrname=p.chr[chr];
		int chrlen;
		if(a.mfn)						// If minus strand is used
		{
			p.getChr(pl,chrname);
			m.getChr(mn,chrname);
			chrlen=pl.size();
			if(chrlen>mn.size()) chrlen=mn.size();
		}
		else								// If minus strand is not used
		{
			p.getChr(pl,chrname);
			mn=pl;
			chrlen=pl.size();
		}
		cout<<"\tLoaded "<<chrname<<", length : "<<chrlen/1000<<"kb"<<endl;
			
		// Read data for gene reads			
		//cout<<"\tReading gene data"<<endl;
		int arraylen=a.mar5+a.mar3;
		for(int i=0;i<g.n;i++)		// Go through all genes
		{
			if(g.chrom[i]==chrname)	// Gene on the same chr
			{
				int s=g.start[i];
				int e=g.end[i];
				if(g.strand[i]=="+")	// plus strand
				{
					for(int j=0;j<arraylen;j++)
					{
						int pos=s+j-a.mar5;
						if(pos>=0&&pos<chrlen)
						{
							pden[i][j]=pl[pos];
							if(mn[pos]>0) mden[i][j]=-mn[pos];
							else mden[i][j]=mn[pos];
						}
					}
				}
				else					// minus strand
				{
					for(int j=0;j<arraylen;j++)
					{
						int pos=e-j+a.mar5;
						if(pos>=0&&pos<chrlen)
						{
							if(mn[pos]>0) pden[i][j]=mn[pos];
							else pden[i][j]=-mn[pos];
							mden[i][j]=-pl[pos];
						}
					}
				}
			}
		}
	}
};

void smooth_density(vector<vector<float> > &den,smoothFunction &sm,int nth,int genegroup)
{
	int n=den.size();
	for(int i=0;i<n;i++)
	{
		if((i%nth)!=genegroup) continue;
		else sm.apply(den[i]);
	}
};
	

