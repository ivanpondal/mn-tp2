#include "mini_test.h"
#include "page_rank.h"
#include "utils.h"

#include <string>
#include <sstream>
#include <math.h>
#include <ctime>
#include <chrono>

void test_cargar_SNAP() {
	string entrada = "tests/test1.txt";
	vector< vector<int> > A = Utils::cargarSNAP(entrada.c_str());
	Utils::imprimirMatriz(A);
}

void test_page_rank_1() {
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
