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
#include <cmath>

using namespace std;

// Parameter initialization
par::par()
{
    mfn= NULL;
    nth= 6;
    // window = 20;
    shStart = -240;
    shEnd = 240;
    bin = 5;
    cutoff =5;
};

// Parsing arguments
bool par::get(int argc, char *argv[])
{
	arg::push("-pl","Plus strand bedgraph filename",pfn);
	arg::push("-mn","Minus bedgraph filename",mfn);
	arg::push("-out","Output file prefix",ofn);
	// arg::push("-win","TSS window size (default=20)",window,true);
	arg::push("-minsh","Minimum divergent shift range (default=-240, negative values for convergent)",shStart,true);
	arg::push("-maxsh","Maximum divergent shift range (default=240)",shEnd,true);
	arg::push("-bin","Bin size (default=5)",bin,true);
	arg::push("-cutoff","Cut-off read count per bin (default=5)",cutoff,true);
	arg::push("-nth","Number of threads (default=6)",nth,true);
	return(arg::get(argc,argv));
};

void get_bgdata_thread(float_2d &pden, float_2d &mden, bgdata &p, bgdata &m, par &a, int chrGroup)
{
	for(int chr=0;chr<p.chr.size();chr++)		// Go through all chromosomes
    {
		if((chr%a.nth)!=chrGroup) continue;     // For multithreading, skip unassigned chromsome
		string chrname=p.chr[chr];
		p.getChr(pden[chr],chrname,a.bin);
		if(m.getChrID(chrname)!=-1)
			m.getChr(mden[chr],chrname,a.bin);
		int chrlen=pden[chr].size();
		if(chrlen<mden[chr].size()) chrlen=mden[chr].size();
        pden[chr].resize(chrlen);
        mden[chr].resize(chrlen);
		for(int i=0;i<chrlen;i++) mden[chr][i]=-mden[chr][i];
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

void make_shbed_thread(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, vector<string> &chromName, par
&a, int chrGroup)
{
	int nchr=pden.size();
    int shBinStart=a.shStart/a.bin;
    int shBinEnd=a.shEnd/a.bin;
    
    for(int chr=0;chr<nchr;++chr)
    {
        if((chr%a.nth)!=chrGroup) continue;
        int chrsize=pden[chr].size();
        int lim=(-shBinStart>shBinEnd)? -shBinStart: shBinEnd;
        for(int i=lim;i<chrsize-lim;++i)
        {
            bedTrack t;
            double maxLgReadSum=0;
            int maxPpos, maxMpos;
            for(int j=shBinStart;j<shBinEnd;++j)
            {
                int ppos,mpos,pc,mc;
                if(j%2)
                {
                    ppos=i+j/2+1;
                    mpos=i-j/2;
                    pc=(int)pden[chr][ppos];
                    mc=(int)mden[chr][mpos];
                }
                else
                {
                    ppos=i+j/2;
                    mpos=i-j/2;
                    pc=(int)pden[chr][ppos];
                    mc=(int)mden[chr][mpos];
                }
                if(pc<a.cutoff||mc<a.cutoff) continue;
                double lgReadSum=log10(pc+1)+log10(mc+1);
                if(lgReadSum>maxLgReadSum)
                {
                    maxLgReadSum=lgReadSum;
                    maxPpos=ppos;
                    maxMpos=mpos;
                }
            }
            if(maxLgReadSum>0) 
            {
                stringstream ss;
                t.chrom=chromName[chr];
                if(maxPpos>maxMpos)
                {
                    t.chromStart=maxMpos*a.bin;
                    t.chromEnd=maxPpos*a.bin;
                    ss<<"div:"<<t.chromEnd-t.chromStart;
                }
                else
                {
                    t.chromStart=maxPpos*a.bin;
                    t.chromEnd=maxMpos*a.bin;
                    ss<<"conv:"<<t.chromEnd-t.chromStart;
                }
                t.name=ss.str();
                t.pair_to_score(pden[chr][maxPpos],mden[chr][maxMpos]);
                if(pden[chr][maxPpos]>mden[chr][maxMpos]) t.strand='+';
                else t.strand='-';
                den[chr].push_back(t);
            }
        }
    }
};

void make_shbed(vector<vector<bedTrack> > &den, float_2d &pden, float_2d &mden, vector<string> &chromName, par &a)
{
	vector<thread> thr(a.nth);
    for(int i=0;i<a.nth;i++)
		thr[i]=thread(make_shbed_thread,ref(den),ref(pden),ref(mden),ref(chromName),ref(a), i);
	for(int i=0;i<a.nth;i++)
		thr[i].join();
};

void save_bed(vector<vector<bedTrack> > &den, string outputPrefix)
{
    stringstream ss;
    ss<<outputPrefix<<".bed";
    ofstream out((char *)ss.str().c_str());

    int nchr=den.size();
    for(int chr=0;chr<nchr;++chr)
        for(int j=0;j<den[chr].size();++j)
            den[chr][j].write(out,6);
};
