#include "mini_test.h"
#include "page_rank.h"
#include "page_rank_esparso.h"
#include "in_deg.h"
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
	vector< vector<int> > A = cargarSNAP(entrada.c_str());
	cout << endl << "Matriz Adyacencia: " << entrada << endl;
	imprimirMatriz(A);
}

void test_matriz_transicion() {
	string entrada = "tests/test1.txt";
	vector< vector<int> > A = cargarSNAP(entrada.c_str());
	cout << endl << "Matriz Transicion: " << entrada << endl;
	imprimirMatriz(matrizTransicion(A));
}

void test_page_rank_1() {
	string entrada = "tests/test1.txt";
	vector< vector<int> > A = cargarSNAP(entrada.c_str());
	PageRank page_rank(matrizTransicion(A));
	vector< PageRank::rankeable > ranking = page_rank.rankear();
	cout << endl << "PageRank: " << entrada << endl;
	for (unsigned int i = 0; i < ranking.size(); i++) {
		cout << ranking[i].posicion << ": " << ranking[i].valor << endl;
	}
}

void test_in_deg_1() {
	string entrada = "tests/test1.txt";
	vector< vector<int> > A = cargarSNAP(entrada.c_str());
	InDeg inDeg(A);
	vector< InDeg::rankeable > ranking = inDeg.rankear();
	cout << endl << "InDeg: " << entrada << endl;
	for (unsigned int i = 0; i < ranking.size(); i++) {
		cout << ranking[i].posicion << ": " << ranking[i].valor << endl;
	}
}

void test_page_rank_esparso_1() {
	string entrada = "tests/test1.txt";
	PageRankEsparso page_rank(cargarSNAPEsparso(entrada.c_str()));
	vector< PageRankEsparso::rankeable > ranking = page_rank.rankear();
	cout << endl << "PageRankEsparso: " << entrada << endl;
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
		test_matriz_transicion();
		test_page_rank_1();
		test_in_deg_1();
		test_page_rank_esparso_1();
	}
	return 0;
}
