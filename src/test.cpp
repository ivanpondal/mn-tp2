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
static int offset = 0;

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
			PageRankEsparso page_rank(cargarSNAPEsparso(path, offset));
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
            // octave --eval "source('footbal_rankings.m'); GeM(exp/gem_resultados_1.txt, exp/gem_resultados_1.out, c=0.85, pres=0.1)"
            char command[1024];
            string out_file(path);
			out_file = remove_extension(out_file);
			out_file.append(".out");
            sprintf(command, "octave --eval \"source('footbal_rankings.m'); GeM('%s','%s', team_codes_filename='data/team_code.txt', c=%.3f, date_limit=0, pres=%.9f);\" >> /dev/null",
                path, out_file.c_str(), tele, tolerancia);
            if(system(command)) { cout << "System failed" << endl; };
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

	vector< map<int, double> > A = cargarSNAPEsparso(in,0);

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

void exp_prank_tiempos_aux(const char * in, FILE * out, double tel, double pres) {
	start_timer();
	resolver(0, tel, 0, in, pres);
	double tiempo = stop_timer();
	fprintf(out, "%s %.2f %.3f %.9f\n", in, tiempo, tel, pres);
}

void exp_prank_tiempos() {
	FILE *file = fopen("exp/pr-tiempos-1-1.out", "w+");
	fprintf(file, "instancia tiempo tel precision\n");

    offset = 0;
    // ~ 6000 nodos
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", file, 0.85, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", file, 0.85, 0.00001);
    exp_prank_tiempos_aux("exp/pr-1-1-p2p-Gnutella08.txt", file, 0.85, 0.0000001);

    // ~ 10000 nodos
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", file, 0.85, 0.001);
	exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", file, 0.85, 0.00001);
    exp_prank_tiempos_aux("exp/pr-1-2-p2p-Gnutella04.txt", file, 0.85, 0.0000001);

    // ~ 36000 nodos
	exp_prank_tiempos_aux("exp/p2p-Gnutella30.txt", file, 0.85, 0.001);
	exp_prank_tiempos_aux("exp/p2p-Gnutella30.txt", file, 0.85, 0.00001);
    exp_prank_tiempos_aux("exp/p2p-Gnutella30.txt", file, 0.85, 0.0000001);

    // ~ 62000 nodos
	exp_prank_tiempos_aux("exp/p2p-Gnutella31.txt", file, 0.85, 0.001);
	exp_prank_tiempos_aux("exp/p2p-Gnutella31.txt", file, 0.85, 0.00001);
    exp_prank_tiempos_aux("exp/p2p-Gnutella31.txt", file, 0.85, 0.0000001);

    offset = 1;
    // ~ 280000 nodos
    exp_prank_tiempos_aux("exp/web-Stanford.txt", file, 0.85, 0.001);
    exp_prank_tiempos_aux("exp/web-Stanford.txt", file, 0.85, 0.00001);
    exp_prank_tiempos_aux("exp/web-Stanford.txt", file, 0.85, 0.0000001);

	fclose(file);
}

void exp_prank_calidad_aux(const char * in, const char * out) {
	FILE *file = fopen(out, "w+");

	vector< vector<int> > A = cargarSNAP(in,0);
	vector< vector<double> > M = matrizTransicion(A);

	fprintf(file, "Matriz Adyacencia\n");
	for (unsigned int i = 0; i < A.size(); ++i) { for (unsigned int j = 0; j < A.size(); ++j) {fprintf(file, "%d ", A[i][j]);} fprintf(file, "\n"); }
	fprintf(file, "\n");

	fprintf(file, "Matriz Transicion\n");
	for (unsigned int i = 0; i < A.size(); ++i) { for (unsigned int j = 0; j < A.size(); ++j) {fprintf(file, "%.6f ", M[i][j]);} fprintf(file, "\n"); }
	fprintf(file, "\n");

	fprintf(file, "Ranking InDeg:\n");
	fprintf(file, "rank nodo valor\n");
	InDeg indeg(A);
	vector< InDeg::rankeable > ranking1 = indeg.rankear();
	for (unsigned int i = 0; i < A.size(); ++i) { fprintf(file, "%d %d %.6f \n", i, ranking1[i].posicion, ranking1[i].valor); }
	fprintf(file, "\n");

	fprintf(file, "Ranking PageRank:\n");
	fprintf(file, "rank nodo valor\n");
	PageRank page_rank(M);
	vector< PageRank::rankeable > ranking2 = page_rank.rankear();
	for (unsigned int i = 0; i < A.size(); ++i) { fprintf(file, "%d %d %.6f \n", i, ranking2[i].posicion, ranking2[i].valor); }
	fprintf(file, "\n");

}

void exp_prank_calidad() {
	set_teletransportacion(0.1);
	exp_prank_calidad_aux("exp/pr-2-1-generated.txt","exp/pr-2-1-1.out");
	set_teletransportacion(0.5);
	exp_prank_calidad_aux("exp/pr-2-1-generated.txt","exp/pr-2-1-2.out");
	set_teletransportacion(0.9);
	exp_prank_calidad_aux("exp/pr-2-1-generated.txt","exp/pr-2-1-3.out");
	exp_prank_calidad_aux("exp/pr-2-2-no-edges.txt","exp/pr-2-2.out");
	exp_prank_calidad_aux("exp/pr-2-3-complete.txt","exp/pr-2-3.out");
}

void exp_gem_resultados() {
    resolver(0, 0, 1, "exp/gem_resultados_1_1.txt", 0.0001);
    resolver(0, 0.3, 1, "exp/gem_resultados_1_2.txt", 0.0001);
    resolver(0, 0.6, 1, "exp/gem_resultados_1_3.txt", 0.0001);
    resolver(0, 0.85, 1, "exp/gem_resultados_1_4.txt", 0.0001);
    resolver(0, 1, 1, "exp/gem_resultados_1_5.txt", 0.0001);

    resolver(0, 0, 1, "exp/gem_resultados_2_1.txt", 0.0001);
    resolver(0, 0.3, 1, "exp/gem_resultados_2_2.txt", 0.0001);
    resolver(0, 0.6, 1, "exp/gem_resultados_2_3.txt", 0.0001);
    resolver(0, 0.85, 1, "exp/gem_resultados_2_4.txt", 0.0001);
    resolver(0, 1, 1, "exp/gem_resultados_2_5.txt", 0.0001);
}

// para correr un test: ./test test.in test.expected {0: EG, 1: LU}
int main(int argc, char *argv[])
{
	srand (time(0));
	// si no hay argumentos corro tests unitarios, si no los de la cátedra
	if(argc == 4){
		int a;
		a = 0;
	}
	else{
		// test_cargar_SNAP();
		// test_matriz_transicion();
		// test_page_rank_1();
		// test_in_deg_1();
		// test_page_rank_esparso_1();

		// exp_prank_manhattan();
		// exp_prank_tiempos();
		// exp_prank_calidad();

        // exp_gem_resultados();
	}
	return 0;
}
