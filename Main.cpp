
#include "EigenvalueSolver.h"


int main(int argc, char const* argv[]) {
	int N, func;
	string filename;
	
	cout << "Enter dimension N of matrix A" << endl;
	cin >> N;
	cout << "Enter problem to solve: 1 for buckling beam, 2 for QM problem - one electron, 3 for QM problem - two electrons " << endl;
	cin >> func;
	//the code spends too long time running when we set N over 100;
	if (N > 100) {
		cout << "N is too large, choose N less than 100" << endl;
		exit(1);
	}
	Test test;
	double conv = 0.0001;
	JacobiSolver jcbs;
	vec eval;
	mat evec;
	double pot =0;
	mat A = zeros<mat>(N, N);
	mat v = eye<mat>(N, N);
	vec r(N) ,ex(N);
	vec wr;
	EigenvalueSolver matrix;
	if (func == 1) {
		filename = "Beam_N_" + to_string(N) + ".txt";
		test.EigenvalueTest();
		matrix.TridiagonalMatrix(A, N);
		double d = A(0, 0);
		double a = A(0, 1);
		test.Exact(N, d, a, ex);
		jcbs.jsolve(N, conv, A, r, v);
		vec eigenvalue = jcbs.get_eigenvals(A, N);
		mat eigenvector = jcbs.get_eigenvecs(A, v, N);
		test.orthogonal(eigenvector);
		eig_sym(eval, evec, A);
		matrix.write_to_file(filename, eigenvalue, ex,eval, N);

	}
	else if(func == 2) {
		matrix.QM_Matrix(A, N, r, func,pot);
		filename = to_string(N) + "N_QM_10.txt";
		test.Exact_Electron(N, ex, func, pot,N);
		jcbs.jsolve(N, conv, A, r, v);
		vec eigenvalue = jcbs.get_eigenvals(A, N);
		mat eigenvector = jcbs.get_eigenvecs(A, v, N);
		test.orthogonal(eigenvector);
		eig_sym(eval, evec, A);
		matrix.write_to_file(filename, eigenvalue, ex, eval, N);
	}
	else if(func == 3) {
		vec ground_values(4);
		wr = { 0.01, 0.5, 1., 5. };
		vec ex2(4);
		vec armad(4);
		for (int i = 0; i < 4; i++) {
			mat A = zeros<mat>(N, N);
			mat v = eye<mat>(N, N);
			pot = wr(i);
			matrix.QM_Matrix(A, N, r, func, pot);
			filename = to_string(N) + "N_QM2_r10.txt";
			test.Exact_Electron(N, ex2, func, pot, i);
			jcbs.jsolve(N, conv, A, r, v);
			vec eigenvalue = jcbs.get_eigenvals(A, N);
			mat eigenvector = jcbs.get_eigenvecs(A, v, N);
			test.orthogonal(eigenvector);
			eig_sym(eval, evec, A);
			armad(i) = eval(0);
		
			ground_values(i) = eigenvalue(0);
			
		}
		matrix.write_to_file_2(filename, N, ex2, wr, ground_values, armad);
		
		}
	else {
		cout << "No valid problem to solve has been chosen" << endl;
		exit(2);
	}


	return 0;
}