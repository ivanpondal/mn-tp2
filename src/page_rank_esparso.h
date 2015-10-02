#ifndef PAGE_RANK_ESPARSO_H_
#define PAGE_RANK_ESPARSO_H_

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

class PageRankEsparso {
	public:
		struct rankeable {
			rankeable(int pos, double v) : posicion(pos), valor(v) {};
			bool operator<(const rankeable& other) const { return valor > other.valor; }

			int posicion;
			double valor;
		};

		PageRankEsparso(vector< map<int, double> > A);
		vector< rankeable > rankear();
		void set_precision(double p);
	private:
		vector< map<int, double> > A;
		double precision = PRECISION_DEFAULT;
		vector<double> metodoPotencia();
		vector<double> multiplyEsparso(const vector< map<int, double> > &A, const vector<double> &x);
};

#endif // PAGE_RANK_ESPARSO_H_INCLUDED
