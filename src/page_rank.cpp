#include <page_rank.h>

PageRank::PageRank(vector< vector<double> > A){
	if (A.size() == 0 || A[0].size() == 0 || A.size() != A[0].size()) {
		throw invalid_argument("PageRank instanciado con una matriz que no es cuadrada");
	}
	this->A = A;
}

void PageRank::set_precision(double p) {
	precision = p;
}

vector< typename PageRank::rankeable > PageRank::rankear() {
	vector<double> autovector = metodoPotencia(); // autovector asociado al autovalor 1

	// cout << "Norma2: " << norma2(autovector) << endl;
	vector< typename PageRank::rankeable > ranking;
	for (unsigned int i = 0; i < autovector.size(); i++) {
		ranking.push_back(rankeable(i, autovector[i]));
	}

	sort(ranking.begin(), ranking.end());
	return ranking;
}

vector<double> PageRank::metodoPotencia() {
	int n = A.size();
	// creo el x aleatorio
	vector<double> x(n, 0);
	for (int i = 0; i < n; i++) {
		x[i] = random_in_range(1,50);
	}

	// traspongo x para poder multiplicarlo por A
	vector< vector<double> > aux = row2Column(x);

	vector< vector<double> > B = multiply(A, A);
	vector< vector<double> > C = A;
	double delta = phi(column2Row(multiply(B, aux))) / phi(column2Row(multiply(C, aux)));
	double last_delta = INFINITY;

	while (fabs(delta - last_delta) > precision) {
		C = B;
		B = multiply(B, A);
		last_delta = delta;
		delta = phi(column2Row(multiply(B, aux))) / phi(column2Row(multiply(C, aux)));
	}
	vector<double> v = column2Row(multiply(B, aux));

	return scaleVector(v, 1/norma1(v));
}

double PageRank::phi(const vector<double> &x) {
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
