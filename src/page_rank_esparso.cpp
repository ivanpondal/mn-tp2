#include <page_rank_esparso.h>

PageRankEsparso::PageRankEsparso(vector< map<int, double> > A){
	this->A = A;
}

void PageRankEsparso::set_precision(double p) {
	precision = p;
}

vector< typename PageRankEsparso::rankeable > PageRankEsparso::rankear() {
	vector<double> autovector = metodoPotencia(); // autovector asociado al autovalor 1

	// cout << "Norma2: " << norma2(autovector) << endl;
	vector< typename PageRankEsparso::rankeable > ranking;
	for (unsigned int i = 0; i < autovector.size(); i++) {
		ranking.push_back(rankeable(i, autovector[i]));
	}

	sort(ranking.begin(), ranking.end());
	return ranking;
}

vector<double> PageRankEsparso::metodoPotencia() {
	int n = A.size();
	// creo el x aleatorio
	vector<double> x(n, 0);
	for (int i = 0; i < n; i++) {
		x[i] = random_in_range(1,50);
	}
	// creo el v aleatorio
	vector<double> v(n, 1/n);

	vector<double> y = x;
	double delta = INFINITY;
	do {
		x = y;
		y = scaleVector(multiplyEsparso(A, x), teletransportacion);
		double w = norma1(x) - norma1(y);
		y = sumVector(y, scaleVector(v, w));

		delta = norma1(sumVector(y, scaleVector(x, double(-1))));
	} while (delta > precision);

	return scaleVector(y, 1/norma2(y));
}

vector<double> PageRankEsparso::multiplyEsparso(const vector< map<int, double> > &A, const vector<double> &x) {
	vector<double> y(x.size());
	for (unsigned int i = 0; i < x.size(); i++) {
		typedef map<int, double>::const_iterator it_type;
		double sum = 0;
		for(it_type iterator = A[i].begin(); iterator != A[i].end(); iterator++) {
			int j = iterator->first;
			sum += (iterator->second)*x[j];
		}
		y[i] = sum;
	}
	return y;
}
