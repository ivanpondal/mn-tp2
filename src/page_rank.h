#ifndef PAGE_RANK_H_
#define PAGE_RANK_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <utils.h>
#include <stdexcept>

using namespace std;
using namespace utils;

#define PRECISION_DEFAULT 0.001

template <class T>
class PageRank{
	public:
		PageRank(vector<vector<T> > A);
		vector< pair<int, T> > rankear();
		void set_precision(double p);
	private:
		vector<vector<T> > A;
		double precision = PRECISION_DEFAULT;
		double metodoPotencia();
		double phi(const vector<T> &A);
};

#endif // PAGE_RANK_H_INCLUDED
