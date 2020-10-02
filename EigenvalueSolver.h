
#ifndef EigenvalueSolver_H
#define EigenvalueSolver_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <time.h>
#include <vector>


using namespace std;
using namespace arma;

class EigenvalueSolver {

public:
	void write_to_file(string filename, vec eigenvalue, vec ex, vec eval, int N);
	void write_to_file_2(string filename, int N, vec ex, vec wr, vec ground_values, vec armad);
	void TridiagonalMatrix(mat&, int N);
	void QM_Matrix(mat&, int N, vec&, int func, double pot);
};
class Test : public EigenvalueSolver {
public:


	void EigenvalueTest();
	void Exact(int N, double d, double a, vec& ex);
	void Exact_Electron(int N, vec& , int func, double pot, int i);
	void orthogonal(mat eigenvec);
};

class JacobiSolver : public EigenvalueSolver {
public:
	
	void jsolve(int,double , mat&, vec&, mat&);
	void find_max(mat , int& , int& , double& , int );
	void print_vals(mat, mat, int, double);
	vector<double> get_eigenvals(mat&, int);
	mat get_eigenvecs(mat&, mat&, int);
	
	
};



#endif