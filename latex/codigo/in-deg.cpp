#include <in_deg.h>

InDeg::InDeg(vector< vector<int> > A){
	if (A.size() == 0 || A[0].size() == 0 || A.size() != A[0].size()) {
		throw invalid_argument("InDeg instanciado con una matriz que no es cuadrada");
	}
	this->A = A;
}

vector< typename InDeg::rankeable > InDeg::rankear() {
	vector< typename InDeg::rankeable > ranking;

	for (unsigned int i = 0; i < A.size(); i++) {
		int grado = 0;
		for (unsigned int j = 0; j < A.size(); j++) {
			if(A[i][j]) {
				grado++;
			}
		}
		ranking.push_back(rankeable(i, grado));
	}

	sort(ranking.begin(), ranking.end());
	return ranking;
}
