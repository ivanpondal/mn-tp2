#ifndef IN_DEG_H_
#define IN_DEG_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <utils.h>
#include <stdexcept>

using namespace std;
using namespace utils;

class InDeg {
	public:
		struct rankeable {
			rankeable(int pos, double v) : posicion(pos), valor(v) {};
			bool operator<(const rankeable& other) const { return valor > other.valor; }

			int posicion;
			double valor;
		};
		InDeg(vector< vector<int> > A);
		vector< rankeable > rankear();
	private:
		vector< vector<int> > A;
};

#endif // IN_DEG_H_INCLUDED
