#include "bg.h"					// Custom reading bedgraph files
#include "arguments.h"			// Parsing command line arguments
#include "smooth.h"				// Custom Gaussian smoothing functions

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

// Declarations
struct args		// Arguments
{
	char *pfn,*mfn;			// Bedgraph filenames
	char *ofn;				// Output filename
	int nth;
	int res,sh,bin;				// resolution
	// Initialization
	args() { nth=6; res=50; sh=80; bin=5; };
	bool get(int argc, char *argv[]);	// Function to parse arguments
};

// Functions used in main
void get_data(vector<vector<float> > &pden,vector<vector<float> > &mden, vector<vector<float> > &den,bgdata &p,bgdata &m,args &a,int chrgroup);
void get_peak(vector<vector<float> > &den,vector<vector<int> > &peakpos, args &a, int chrgroup);

int main(int argc,char *argv[])
{
	// Parsing arguments
	args a;
	if(!a.get(argc,argv)) return 0;
	
	// Load bedgraph data
	bgdata p,m;
	p.load(a.pfn);		// Plus strand data
	m.load(a.mfn);	// Minus strand data if provided
	int nchr=p.chr.size();
		
	// Arrays of bedgraph data values for every gene
	vector<vector<float> > pden(nchr),mden(nchr),den(nchr);	// Data arrays for bedgraph values
	// Read data into data array, retrieve read counts
	// Use multithreading
	vector<thread> thr(a.nth);
	cout<<"Reading bedgraphs"<<endl;
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_data,ref(pden),ref(mden),ref(den),ref(p),ref(m),ref(a),i);
	for(int i=0;i<a.nth;i++)
	{
		thr[i].join();
	}
	cout<<"Finished reading bedgraphs"<<endl;

	vector<vector<int> > ppeakpos(nchr);
	vector<vector<int> > mpeakpos(nchr);
	vector<vector<int> > peakpos(nchr);
	cout<<"Scanning peaks"<<endl;
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_peak,ref(pden),ref(ppeakpos),ref(a),i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_peak,ref(mden),ref(mpeakpos),ref(a),i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_peak,ref(den),ref(peakpos),ref(a),i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
	cout<<"Finished scanning peaks"<<endl;

	ofstream out(a.ofn);
	out<<"chr\tpos\tminus_pos\tplus_pos\tminus_val\tplus_val"<<endl;
	int sb=a.sh/a.bin;
	for(int i=0;i<peakpos.size();i++)
	{
		if(ppeakpos[i].size()==0||mpeakpos[i].size()==0) continue;
		for(int j=0;j<peakpos[i].size();j++)
		{
			
			int pi,p1,p2,mi,m1,m2,pp=peakpos[i][j];
			pi=upper_bound(ppeakpos[i].begin(),ppeakpos[i].end(),pp+sb)-ppeakpos[i].begin();
			if(pi==ppeakpos[i].size()) p1=p2=ppeakpos[i][pi-1];
			else if(pi==0) p1=p2=ppeakpos[i][0];
			else
			{
				p1=ppeakpos[i][pi];
				p2=ppeakpos[i][pi-1];
			}
			mi=upper_bound(mpeakpos[i].begin(),mpeakpos[i].end(),pp-sb)-mpeakpos[i].begin();
			if(mi==mpeakpos[i].size()) m1=m2=mpeakpos[i][mi-1];
			else if(mi==0) m1=m2=mpeakpos[i][0];
			else
			{
				m1=mpeakpos[i][mi];
				m2=mpeakpos[i][mi-1];
			}
			if(p2<pp) p2=p1;
			if(m1>pp) m1=m2;
			if(p1-pp-sb>pp+sb-p2) p1=p2;
			if(m1-pp+sb>pp-sb-m2) m1=m2;
			out<<p.chr[i]<<"\t"<<peakpos[i][j]*a.bin<<"\t"<<m1*a.bin<<"\t"<<p1*a.bin<<"\t";
			out<<mden[i][m1]<<"\t"<<pden[i][p1]<<endl;
		}
	}
}
	   
// Parsing arguments
bool args::get(int argc, char *argv[])
{
	// Parsing arguments
	arg::push("-p","plus bedgraph filename",pfn);
	arg::push("-m","minus bedgraph filename (optional)",mfn,true);
	arg::push("-o","output filename",ofn);
	arg::push("-r","peak resolution (default=50)",res,true);
	arg::push("-b","bin size (default=5)",bin,true);
	arg::push("-sh","shift (default=80)",sh,true);
	arg::push("-nth","Number of threads (default=6)",nth,true);
	return(arg::get(argc,argv));
};

// Function for making a data array from begraph files
void get_data(vector<vector<float> > &pden,vector<vector<float> > &mden,vector<vector<float> > &den,bgdata &p,bgdata &m,args &a,int chrgroup)
{
	double sw=a.res/a.bin;
	smoothFunction sf(sw);
	int sb=a.sh/a.bin;
	for(int chr=0;chr<p.chr.size();chr++)		// Go through all chromosomes
	{
		if((chr%a.nth)!=chrgroup) continue;	// For multithreading
		string chrname=p.chr[chr];
		p.getlevel(pden[chr],chrname,0,p.end[p.getChrID(chrname)].back(),a.bin);
		if(m.getChrID(chrname)!=-1)
			m.getlevel(mden[chr],chrname,0,m.end[m.getChrID(chrname)].back(),a.bin);
		int chrlen=pden[chr].size();
		if(chrlen<mden[chr].size()) chrlen=mden[chr].size();
		if(chrlen<sw) chrlen=sw;
		pden[chr].resize(chrlen);
		mden[chr].resize(chrlen);
		for(int i=0;i<chrlen;i++) mden[chr][i]=-mden[chr][i];
		sf.apply(pden[chr]);
		sf.apply(mden[chr]);
		den[chr].resize(chrlen);
		for(int i=sb;i<chrlen-sb;i++)
			if(pden[chr][i+sb]>0&&mden[chr][i-sb]>0)
				den[chr][i]=pden[chr][i+sb]*mden[chr][i-sb];
	}
};

void get_peak(vector<vector<float> > &den,vector<vector<int> > &peakpos, args &a,int chrgroup)
{
	int nchr=den.size();
	for(int chr=0;chr<nchr;++chr)
	{
		if((chr%a.nth)!=chrgroup) continue;
		int chrsize=den[chr].size();
		for(int i=2;i<chrsize-2;++i)
		{
			float x1=den[chr][i-2],x2=den[chr][i-1],x3=den[chr][i],x4=den[chr][i+1],x5=den[chr][i+2];
			if(x3>=x2&&x3>=x1&&x3>=x4&&x3>=x5&&!(x3==x1&&x3==x5))
				peakpos[chr].push_back(i);
		}
	}
};
