#ifndef _BEDGRAPH_H
#define _BEDGRAPH_H

#include <string>
#include <vector>

struct bgdata
{
    int load(char *fn);
    void getlevel(std::vector<float> &o, std::string &c, int s, int e, int r=1, bool isaverage=false);

    std::vector<std::string> chr;
	std::vector<std::vector<int> > start;
    std::vector<std::vector<int> > end;
	std::vector<std::vector<float> > val;
    std::vector<float> v;

    int ch;
    int res;
    int len;
	bool set(std::string &c, int r=1);
    void getcoverage(std::vector<float> &o, std::string &c, int s, int e, int r=1);
    int getChrID(std::string &c);
    void getChr(std::vector<float> &o, std::string &c);
    float &operator[](int i);
    float operator()(int c,int i);
    float operator()(std::string &c,int i);
};
#endif
