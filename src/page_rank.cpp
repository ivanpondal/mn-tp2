#include <page_rank.h>

template <class T>
PageRank<T>::PageRank(vector<vector<T> > A){
	if (A.size() == 0 || A[0].size() == 0 || A.size() != A[0].size()) {
		throw invalid_argument("PageRank instanciado con una matriz que no es cuadrada");
	}
	this->A = A;
}

template <class T>
void PageRank<T>::set_precision(double p) {
	precision = p;
}

template <class T>
vector< pair<int, T> > PageRank<T>::rankear() {
	double autovalor = metodoPotencia();
	cout << autovalor << endl;
	return vector< pair<int, T> >();
}

template <class T>
double PageRank<T>::metodoPotencia() {
	int n = A.size();
	// creo el x aleatorio
	vector<T> x(n, 0);
	for (int i = 0; i < n; i++) {
		x[i] = random_in_range(1,50);
	}

	// traspongo x para poder multiplicarlo por A
	vector< vector<T> > aux = row2Column(x);

	vector< vector<T> > B = multiply(A, A);
	vector< vector<T> > C = A;
	double delta = phi(column2Row(multiply(B, aux))) / phi(column2Row(multiply(C, aux)));
	double last_delta = INFINITY;

	while (fabs(delta - last_delta) > precision) {
		C = B;
		B = multiply(B, A);
		last_delta = delta;
		delta = phi(column2Row(multiply(B, aux))) / phi(column2Row(multiply(C, aux)));
	}

	return delta;
}

template <class T>
double PageRank<T>::phi(const vector<T> &x) {
	T ret = x[0];
	T max = abs(ret);
	for (int i = 1; i < x.size(); i++) {
		if (max < abs(x[i])) {
			ret = x[i];
			max = abs(ret);
		}
	}

	return ret;
}

// pongo explicito cuales son las posibles instancias de page_rank:
template class PageRank<double>;
template class PageRank<int>;
