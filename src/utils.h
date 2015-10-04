#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <iomanip>
#include <stdlib.h>
#include <ctime>

using namespace std;

#define TELETRANSPORTACION_DEFAULT 0.85

namespace utils {
	static double teletransportacion = TELETRANSPORTACION_DEFAULT;

	static int random_in_range(int min, int max) {
		srand (time(NULL));
		return min + (rand() % (max - min + 1));
	}


	static void set_teletransportacion(double t) {
		teletransportacion = t;
	}

	template <typename T>
	static vector< vector<T> > vector2matrix(const vector<T> &x) {
		vector< vector<T> > aux;
		aux.push_back(x);
		return aux;
	}

	template <typename T>
	static vector< vector<T> > row2Column(const vector<T> &x) {
		vector< vector<T> > aux(x.size());
		for (unsigned int i = 0; i < x.size(); ++i) {
			vector<T> elem;
			elem.push_back(x[i]);
			aux[i] = (elem);
		}
		return aux;
	}

	template <typename T>
	static vector<T> column2Row(const vector< vector<T> > &x) {
		vector<T> aux(x.size());
		for (unsigned int i = 0; i < x.size(); ++i) {
			aux[i] = x[i][0];
		}
		return aux;
	}

	static double norma2(const vector<double> &x) {
		double sum = 0;
		for (unsigned int i = 0; i < x.size(); i++) {
			sum += pow(x[i],2);
		}

		return sqrt(sum);
	}

	static double norma1(const vector<double> &x) {
		double sum = 0;
		for (unsigned int i = 0; i < x.size(); i++) {
			sum += fabs(x[i]);
		}
		return sum;
	}

	static double normaManhattan(const vector<double> &x, const vector<double> &y) {
		double sum = 0;
		for (unsigned int i = 0; i < x.size(); i++) {
			sum += fabs(x[i]-y[i]);
		}
		return sum;
	}

	template <typename T>
	static vector< vector<T> > sum(const vector< vector<T> > &A, const vector< vector<T> > &B) {
		if (A.size() !=  B.size() || A[0].size() != B[0].size()) {
			cout << "No coinciden los tamaños para sumar las matrices" << endl;
			cout << "A: " << A.size() << "x" << A[0].size() << endl;
			cout << "B: " << B.size() << "x" << B[0].size() << endl;
			return vector< vector<T> >();
		}

		int rows = A.size();
		int cols = A[0].size();
		vector< vector<T> > resultado(rows, vector<T>(cols, 0));

		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				resultado[i][j] += A[i][j] + B[i][j];
			}
		}
		return resultado;
	}

	template <typename T>
	static vector<T> sumVector(const vector<T> &x, const vector<T> &y) {
		if (x.size() !=  y.size()) {
			cout << "No coinciden los tamaños para sumar los vectores" << endl;
			cout << "x: " << x.size() << endl;
			cout << "y: " << y.size() << endl;
			return vector<T>();
		}

		vector<T> resultado(x.size(), 0);

		for (int i = 0; i < x.size(); ++i) {
			resultado[i] = x[i] + y[i];
		}
		return resultado;
	}

	template <typename T>
	static vector<T> scaleVector(const vector<T> &A, const T &scalar) {
		int n = A.size();
		vector<T> resultado(n, 0);

		for (int i = 0; i < n; ++i) {
			resultado[i] = A[i] * scalar;
		}
		return resultado;
	}

	template <typename T>
	static vector< vector<T> > scaleMatriz(const vector< vector<T> > &A, const T &scalar) {
		int rows = A.size();
		int cols = A[0].size();
		vector< vector<T> > resultado(rows, vector<T>(cols, 0));

		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				resultado[i][j] = A[i][j] * scalar;
			}
		}
		return resultado;
	}

	template <typename T>
	static vector< vector<T> > multiply(const vector< vector<T> > &A, const vector< vector<T> > &B) {
		if (A.size() == 0 || B.size() == 0  || A[0].size() != B.size()) {
			cout << "No coinciden los tamaños para multiplicar las matrices" << endl;
			cout << "A: " << A.size() << "x" << A[0].size() << endl;
			cout << "B: " << B.size() << "x" << B[0].size() << endl;
			return vector< vector<T> >();
		}

		int rows = A.size();
		int cols = B[0].size();
		int dim = B.size();
		vector< vector<T> > resultado(rows, vector<T>(cols, 0));

		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				for (int k = 0; k < dim; ++k) {
					resultado[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		return resultado;
	}

	/**
	* Toma una matriz de conectividad (w_{i,j} tiene un 1 si j tiene un link saliente a i)
	* y devuelve su matriz de transición asociada.
	*/
	static vector< vector<double> > matrizTransicion(const vector< vector<int> > &A) {
		int n = A.size();
		// convierto la matriz A en una matriz P, donde p_{i,j} = 1/n_{j} si w_{i,j} = 1, y 0 en caso contrario.
		// Y n_{j} = cantidad de links salientes desde la página j
		vector< vector<double> >P(n, vector<double>(n, 0));
		vector<int> links_salientes(n,0);
		for (int j = 0; j < n; j++) {
			// cuento los links salientes del nodo j
			for (int i = 0; i < n; i++) {
				links_salientes[j] += A[i][j];
			}

			// asigno el peso a la celda p_{i,j}
			for (int i = 0; i < n; i++) {
				P[i][j] = links_salientes[j] > 0 ? double(A[i][j])/double(links_salientes[j]) : double(1)/double(n);
			}
		}

		// convierto la matriz P en una matriz P1 = P + D, con D = v*d^{t}, con v_{i} = 1/n y d_{i} = 1 si n_{i} = 0, y 0 en caso contrario
		vector<double> v(n, double(1)/double(n));
		vector<double> d(n, 0);
		for (int j = 0; j < n; j++) {
			if (links_salientes[j] == 0) {
				d[j] = 1;
			}
		}

		vector< vector<double> > D = multiply(row2Column(v), vector2matrix(d));
		vector< vector<double> > P1 = sum(P, D);

		// agrego la probabilidad del navegante aleatorio: P2 = c*P1 + (1-c)*E, con E = v*1^{t}
		vector<double> unos(n, 1);
		vector< vector<double> > E = multiply(row2Column(v), vector2matrix(unos));

		vector< vector<double> > P2 = sum(scaleMatriz(P1, teletransportacion), scaleMatriz(E, 1-teletransportacion));

		return P2;
	}

	static vector< vector<int> > cargarSNAP(const char* entrada, int offset = 1) {
		ifstream datos;
		datos.open(entrada);
		if (!datos.good()) {
			cout << endl;
			cout << "\tNo se pudo cargar el archivo: " << entrada << endl;
		}

		// skip first two lines
		string line;
		getline(datos,line);
		getline(datos,line);
		// get how many nodes and edges has the graph
		int nodes = 0, edges = 0;
		datos >> line >> line; // Skip "# Nodes:"
		datos >> nodes;
		datos >> line; // Skip " Edges:"
		datos >> edges;
		// skip another line
		datos  >> line;
		getline(datos,line);
		// load the graph edges
		vector< vector<int> >A(nodes, vector<int>(nodes, 0));
		for (int e = 0; e < edges; e++) {
			int i, j;
			datos >> i >> j;
			A[j-offset][i-offset] = 1;
		}
		datos.close();

		return A;
	}

	static vector< map<int, double> > cargarSNAPEsparso(const char* entrada, int offset = 1) {
		ifstream datos;
		datos.open(entrada);
		if (!datos.good()) {
			cout << endl;
			cout << "\tNo se pudo cargar el archivo: " << entrada << endl;
		}

		// skip first two lines
		string line;
		getline(datos,line);
		getline(datos,line);
		// get how many nodes and edges has the graph
		int nodes = 0, edges = 0;
		datos >> line >> line; // Skip "# Nodes:"
		datos >> nodes;
		datos >> line; // Skip " Edges:"
		datos >> edges;
		// skip another line
		datos  >> line;
		getline(datos,line);
		// load the graph edges
		vector< map<int, double> > A(nodes);
		vector<int> links_salientes(nodes);
		for (int e = 0; e < edges; e++) {
			int i, j;
			datos >> i >> j;
			A[j-offset].insert(pair<int,double>(i-offset,1));
			links_salientes[i-offset]++;
		}
		datos.close();

		// convierto la matriz A en una matriz P, donde p_{i,j} = 1/n_{j} si w_{i,j} = 1, y 0 en caso contrario.
		// Y n_{j} = cantidad de links salientes desde la página j
		for (int i = 0; i < nodes; i++) {
			typedef map<int, double>::iterator it_type;
			for(it_type iterator = A[i].begin(); iterator != A[i].end(); iterator++) {
				int j = iterator->first;
				iterator->second = double(1)/double(links_salientes[j]);
			}
		}

		return A;
	}

	static vector< vector<double> > cargarLigaDeportiva(const char* entrada) {
		return vector< vector<double> >();
	}

	template <typename T>
	static void imprimirMatriz(const vector< vector<T> > &A) {
		for(unsigned int i = 0; i < A.size(); i++){
			cout<<"[";
			for(unsigned int j = 0; j < A[i].size(); j++){
				cout<<" "<< A[i][j];
			}
			cout<<" ]"<<endl;
		}
	}

	template <typename T>
	static void imprimirVector(const vector<T> &A) {
		cout<<"[";
		for(unsigned int i = 0; i < A.size(); i++){
			cout<<" "<< A[i];
		}
		cout<<" ]"<<endl;
	}
};

#endif // UTILS_H_INCLUDED
