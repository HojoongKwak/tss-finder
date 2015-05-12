#include "bed.h"
#include <iostream>

using namespace std;

int main()
{
    bedTrack t;
    while(true)
    {
        cout<<"score?";
        cin>>t.score;
        int a,b;
        t.score_to_pair(a,b);
        cout<<a<<"\t"<<b<<endl;
    }
}

