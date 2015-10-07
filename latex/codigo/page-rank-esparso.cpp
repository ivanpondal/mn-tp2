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
	vector<double> v(n, double(1)/double(n));

	vector<double> y = x;
	double delta = INFINITY;
	double last_delta;
	do {
		x = y;
		y = scaleVector(multiplyEsparso(A, x), teletransportacion);
		double w = norma1(x) - norma1(y);
		y = sumVector(y, scaleVector(v, w));
		last_delta = delta;
		delta = phi(y) / phi(x);
	} while (fabs(delta - last_delta) > precision);

	return scaleVector(y, 1/norma1(y));
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

void PageRankEsparso::imprimirEsparso(const vector< map<int, double> > &A) {
	vector< vector<double> > v(A.size(), vector<double>(A.size(), 0));
	for (unsigned int i = 0; i < A.size(); i++) {
		typedef map<int, double>::const_iterator it_type;
		for(it_type iterator = A[i].begin(); iterator != A[i].end(); iterator++) {
			v[i][iterator->first] = iterator->second;
		}
	}
	imprimirMatriz(v);
}

double PageRankEsparso::phi(const vector<double> &x) {
	double ret = x[0];
	double max = abs(ret);
	for (unsigned int i = 1; i < x.size(); i++) {
		if (max < abs(x[i])) {
			ret = x[i];
			max = abs(ret);
		}
	}

	return ret;
}
