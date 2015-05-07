#include "bg.h"					// Custom reading bedgraph files
#include "smooth.h"				// Custom Gaussian smoothing functions
#include "arguments.h"

#include "find_tss.h"

#include <string>
#include <vector>
#include <thread>
#include <functional>

using namespace std;

// Parameter initialization
par::par()
{
    mfn= NULL;
    nth= 6;
    window = 20;
    shStart = -120;
    shEnd = 120;
    bin =5;
};

// Parsing arguments
bool par::get(int argc, char *argv[])
{
	arg::push("-p","Plus strand bedgraph filename",pfn);
	arg::push("-m","Minus bedgraph filename (optional)",mfn,true);
	arg::push("-o","Output file prefix",ofn);
	arg::push("-w","TSS window size (default=20)",window,true);
	arg::push("-ss","Maximum convergent shift range (default=-120)",shStart,true);
	arg::push("-se","Maximum divergent shift range (default=120)",shEnd,true);
	arg::push("-b","Bin size (default=5)",bin,true);
	arg::push("-nth","Number of threads (default=6)",nth,true);
	return(arg::get(argc,argv));
};

void get_bgdata_thread(float_2d &pden, float_2d &mden, bgdata &p, bgdata &m, par &a, int chrGroup)
{
    int wc=a.window/a.bin;
	for(int chr=0;chr<p.chr.size();chr++)		// Go through all chromosomes
    {
		if((chr%a.nth)!=chrGroup) continue;     // For multithreading, skip unassigned chromsome
        vector<float> praw, mraw;
		string chrname=p.chr[chr];
		p.getChr(praw,chrname,a.bin);
		if(m.getChrID(chrname)!=-1)
			m.getChr(mraw,chrname,a.bin);
		int chrlen=praw.size();
		if(chrlen<mraw.size()) chrlen=mraw.size();
		praw.resize(chrlen);
		mraw.resize(chrlen);
        pden[chr].resize(chrlen);
        mden[chr].resize(chrlen);
		for(int i=0;i<chrlen;i++)
        {
            for(int j=0;j<wc;++j)
            {
                int pos=i-wc/2+j;
                if(pos<0) pos=0;
                else if(pos>=chrlen) pos=chrlen-1;
                pden[chr][i]+=praw[pos];
                mden[chr][i]-=mraw[pos];
            }
        }
	}
};

void get_bgdata(float_2d &pden, float_2d &mden, bgdata &p, bgdata &m, par &a)
{
	vector<thread> thr(a.nth);
	for(int i=0;i<a.nth;i++)
		thr[i]=thread(get_bgdata_thread,ref(pden),ref(mden),ref(p),ref(m),ref(a),i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
};

void make_shbed_thread(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, par &a, int shift, int chrGroup)
{
	int nchr=pden.size();
    int sw=shift/a.bin;
    
    for(int chr=0;chr<nchr;++chr)
    {
        if((chr%a.nth)!=chrGroup) continue;
        int chrsize=pden[chr].size();
        for(int i=0;i<chrsize;++i)
        {
    
        }
    }
};

void make_shbed(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, par &a, int shift)
{
	vector<thread> thr(a.nth);
    for(int i=0;i<a.nth;i++)
		thr[i]=thread(make_shbed_thread,ref(den),ref(pden),ref(mden),ref(a),shift, i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
};

