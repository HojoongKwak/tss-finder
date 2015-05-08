#include "bg.h"					// Custom reading bedgraph files
#include "smooth.h"				// Custom Gaussian smoothing functions
#include "arguments.h"

#include "find_tss.h"

#include <string>
#include <vector>
#include <thread>
#include <sstream>
#include <iostream>
#include <fstream>
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
	arg::push("-pl","Plus strand bedgraph filename",pfn);
	arg::push("-mn","Minus bedgraph filename",mfn);
	arg::push("-out","Output file prefix",ofn);
	arg::push("-win","TSS window size (default=20)",window,true);
	arg::push("-minsh","Minimum divergent shift range (default=-120, negative values for convergent)",shStart,true);
	arg::push("-maxsh","Maximum divergent shift range (default=120)",shEnd,true);
	arg::push("-bin","Bin size (default=5)",bin,true);
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

void make_shbed_thread(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, vector<string> &chromName, par &a, int shift, int chrGroup)
{
	int nchr=pden.size();
    int sw=shift/a.bin;
    stringstream ss;
    ss<<"tsp."<<shift<<".bs";
    string trackName=ss.str();
    for(int chr=0;chr<nchr;++chr)
    {
        if((chr%a.nth)!=chrGroup) continue;
        int chrsize=pden[chr].size();
        for(int i=sw;i<chrsize-sw;++i)
        {
            bedTrack t;
            int pc=(int)pden[chr][i+sw];
            int mc=(int)mden[chr][i-sw];
            if(pc>0&&mc>0)
            {
                t.chrom=chromName[chr];
                t.chromStart=(i-sw)*a.bin;
                t.chromEnd=(i+sw)*a.bin;
                t.name=trackName;
                t.pair_to_score(pc,mc);
                if(pc>mc) t.strand='+';
                else t.strand='-';
                den[chr].push_back(t);
            }
        }
    }
};

void make_shbed(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, vector<string> &chromName, par &a, int shift)
{
	vector<thread> thr(a.nth);
    for(int i=0;i<a.nth;i++)
		thr[i]=thread(make_shbed_thread,ref(den),ref(pden),ref(mden),ref(chromName),ref(a),shift, i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
};

void save_bed(vector<vector<bedTrack> > &den, string outputPrefix, int shift)
{
    stringstream ss;
    ss<<outputPrefix<<"."<<shift<<"bp.bed";
    ofstream out((char *)ss.str().c_str());

    int nchr=den.size();
    for(int chr=0;chr<nchr;++chr)
        for(int j=0;j<den[chr].size();++j)
            den[chr][j].write(out,6);
};
