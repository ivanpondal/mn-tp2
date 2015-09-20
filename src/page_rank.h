#ifndef PAGE_RANK_H_
#define PAGE_RANK_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <utils.h>

using namespace std;

class PageRank{
	public:
		PageRank(vector<vector<double> > A);
		vector< pair<int, double> > rankear();
	private:
		vector<vector<double> > A;
		double metodoPotencia();
};

#endif // PAGE_RANK_H_INCLUDED
