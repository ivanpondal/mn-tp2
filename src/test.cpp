#include "mini_test.h"
#include "page_rank.h"
#include "utils.h"

#include <string>
#include <sstream>
#include <math.h>
#include <ctime>
#include <chrono>

using namespace std;
using namespace utils;

void test_cargar_SNAP() {
	string entrada = "tests/test1.txt";
	vector< vector<double> > A = cargarSNAP(entrada.c_str());
	imprimirMatriz(A);
}

void test_page_rank_1() {
	string entrada = "tests/test1.txt";
	vector< vector<double> > A = cargarSNAP(entrada.c_str());
	PageRank page_rank(A);
	vector< PageRank::rankeable > ranking = page_rank.rankear();
	for (unsigned int i = 0; i < ranking.size(); i++) {
		cout << ranking[i].posicion << ": " << ranking[i].valor << endl;
	}
}

// para correr un test: ./test test.in test.expected {0: EG, 1: LU}
int main(int argc, char *argv[])
{
	// si no hay argumentos corro tests unitarios, si no los de la cátedra
	if(argc == 4){
		int a;
		a = 0;
	}
	else{
		test_cargar_SNAP();
		test_page_rank_1();
	}
	return 0;
}
