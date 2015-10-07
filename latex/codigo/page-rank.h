#ifndef PAGE_RANK_H_
#define PAGE_RANK_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <utils.h>
#include <stdexcept>

using namespace std;
using namespace utils;

#define PRECISION_DEFAULT 0.001

class PageRank{
	public:
		struct rankeable {
			rankeable(int pos, double v) : posicion(pos), valor(v) {};
			bool operator<(const rankeable& other) const { return valor > other.valor; }

			int posicion;
			double valor;
		};

		PageRank(vector< vector<double> > A);
		vector< rankeable > rankear();
		void set_precision(double p);
	private:
		vector< vector<double> > A;
		double precision = PRECISION_DEFAULT;
		vector<double> metodoPotencia();
		double phi(const vector<double> &A);
};

#endif // PAGE_RANK_H_INCLUDED
