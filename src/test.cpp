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

static chrono::time_point<chrono::high_resolution_clock> start_time;

void start_timer() {
    start_time = chrono::high_resolution_clock::now();
}

double stop_timer() {
    chrono::time_point<chrono::high_resolution_clock> end_time = chrono::high_resolution_clock::now();
    return double(chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count());
}

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

// Programa principal que pide la cátedra

string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

void resolver(int algoritmo, double tele, int tipo_instancia, const char* path, double tolerancia) {
	if (algoritmo == 0) {
		// Resolver con PageRank
		if (tipo_instancia == 0) {
			// la instancia es de paginas web
			PageRankEsparso page_rank(cargarSNAPEsparso(path, 0));
			page_rank.set_precision(tolerancia);
			set_teletransportacion(tele);
			vector< PageRankEsparso::rankeable > ranking = page_rank.rankear();
			// imprimo el resultado
			string out_file(path);
			out_file = remove_extension(out_file);
			out_file.append(".out");

			FILE *file = fopen(out_file.c_str(), "w+");
			for (unsigned int i = 0; i < ranking.size(); i++) {
				fprintf(file, "%d %5.f\n",ranking[i].posicion, ranking[i].valor);
			}
		} else {
			// la instancia es de deportes
		}
	} else {
		// Resolver con Método alternativo
	}
}

// Experimentos

double phi(const vector<double> &x) {
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

vector<double> multiplyEsparso(const vector< map<int, double> > &A, const vector<double> &x) {
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

void exp_prank_manhattan_aux(const char * in, const char * out, double tel) {

	double tel_aux = teletransportacion;
	teletransportacion = tel;

	FILE *file = fopen(out, "w+");
	fprintf(file, "instancia l1\n");

	vector< map<int, double> > A = cargarSNAPEsparso(in);

	int n = A.size();

	vector<double> x(n, 0);
	for (int i = 0; i < n; i++) {
		x[i] = random_in_range(1,50);
	}

	// creo el v aleatorio
	vector<double> v(n, double(1)/double(n));

	vector<double> y = x;
	double delta = INFINITY;
	double last_delta;

	int i = 0;
	do {
		x = y;
		y = scaleVector(multiplyEsparso(A, x), teletransportacion);
		double w = norma1(x) - norma1(y);
		y = sumVector(y, scaleVector(v, w));
		last_delta = delta;
		delta = phi(y) / phi(x);

		double manh = normaManhattan(y,x);
		fprintf(file, "%d %.6f\n", i, manh);
		i++;

	} while (fabs(delta - last_delta) > 0.0001);

	fclose(file);

	teletransportacion = tel_aux;
}

void exp_prank_manhattan() {
	exp_prank_manhattan_aux("exp/pr-1-1-p2p-Gnutella08.txt", "exp/pr-1-1-1.out", 0.3);
	exp_prank_manhattan_aux("exp/pr-1-1-p2p-Gnutella08.txt", "exp/pr-1-1-2.out", 0.6);
	exp_prank_manhattan_aux("exp/pr-1-1-p2p-Gnutella08.txt", "exp/pr-1-1-3.out", 0.9);

	exp_prank_manhattan_aux("exp/pr-1-2-p2p-Gnutella04.txt", "exp/pr-1-2-1.out", 0.3);
	exp_prank_manhattan_aux("exp/pr-1-2-p2p-Gnutella04.txt", "exp/pr-1-2-2.out", 0.6);
	exp_prank_manhattan_aux("exp/pr-1-2-p2p-Gnutella04.txt", "exp/pr-1-2-3.out", 0.9);
}

void exp_prank_tiempos_aux(const char * in, double tel, double pres) {
	start_timer();
	resolver(0, tel, 0, in, pres);
	double tiempo = stop_timer();
	cout  << "Tiempo para: "<< in << " (precision: "  << pres << ") -> " << tiempo << " ns" << endl;
}

void exp_prank_tiempos() {
	cout << endl << "Experimentacion tiempos PageRank Esparso" << endl;
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.3, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.6, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.85, 0.001);

	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.3, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.6, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.85, 0.001);

	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.3, 0.00001);
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.6, 0.00001);
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", 0.85, 0.00001);

	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.3, 0.00001);
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.6, 0.00001);
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", 0.85, 0.00001);
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
		// exp_prank_manhattan();
		exp_prank_tiempos();
	}
	return 0;
}
