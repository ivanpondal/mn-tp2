#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <iomanip>

using namespace std;

class Utils {
	public:
		static vector<vector<double> > multiply(const vector<vector<double> > &A, const vector<vector<double> > &B) {
			if (A.size() == 0 || B.size() == 0  || A[0].size() != B.size()) {
				cout << "No coinciden los tamaños para multiplicar las matrices" << endl;
				cout << "A: " << A.size() << "x" << A[0].size() << endl; // 900 x 900
				cout << "B: " << B.size() << "x" << B[0].size() << endl; // 1 x 900
				return vector<vector<double> >();
			}

			int rows = A.size();
			int cols = B[0].size();
			int dim = B.size();
			vector<vector<double> > resultado(rows, vector<double>(cols, 0));

			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					for (int k = 0; k < dim; ++k) {
						resultado[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			return resultado;
		}

		static vector<vector<int> > cargarSNAP(const char* entrada) {
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
			// loead the graph edges
			vector<vector<int> >A(nodes, vector<int>(nodes, 0));
			for (int e = 0; e < edges; e++) {
				int i, j;
				datos >> i >> j;
				A[i-1][j-1] = 1;
			}
			datos.close();

			return A;
		}

		static vector<vector<double> > cargarLigaDeportiva(const char* entrada) {
			return vector<vector<double> >();
		}

		template <typename T>
		static void imprimirMatriz(const vector<vector<T> > &A) {
			for(unsigned int i = 0; i < A.size(); i++){
				cout<<"[";
				for(unsigned int j = 0; j < A[i].size(); j++){
					cout<<" "<< A[i][j];
				}
				cout<<" ]"<<endl;
			}
		}
};

#endif // UTILS_H_INCLUDED
