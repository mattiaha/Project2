#include "EigenvalueSolver.h"



void Test::EigenvalueTest() {
	int N = 3;
	mat B = zeros<mat>(N, N);
	mat v = eye<mat>(N, N);
	vec r(N), ex(N);
	EigenvalueSolver matrix;
	matrix.TridiagonalMatrix(B, N);

	double d = B(0, 0);
	double a = B(0, 1);
	Exact(N, d, a, ex);
	double h = 1. / (N + 1);
	double hh = h * h;
	JacobiSolver jcbs;
	double conv = 0.01;
	jcbs.jsolve(N, conv, B, r, v);
	double pi = 3.141592;
	vec eigen = jcbs.get_eigenvals(B, N);
	double tol = 1E-3;

	for (int i = 0; i < N; i++) {

		if (abs(ex(i) - eigen(i)) > tol) {
			cout << "The generated eigenvalues differ too much from analytical eigenvectors (tolerance is" << tol << ")" << endl;
			exit(3);

		}
	}
}
void Test::Exact(int N, double d, double a, vec& ex) {
	double pi = 3.141592;
	for (int i = 0; i < N; i++) {
		// assume we have a tridiagonal matrix ;
		ex(i) = d + 2 * a * cos(pi * (i + 1) / (N + 1));
	}


}

void Test::Exact_Electron(int N, vec& ex, int func, double pot, int i) {
	if (func == 2) {
		ex(0) = 3.;
		for (int i = 1; i < N; i++) {
			ex(i) = ex(i - 1) + 4.;
		}
	}
	if (func == 3) {
		
		if (pot <= 1.0) {
			
			double V = (3. / 2.) * pow((pot / 2.), (2./3.));
			ex(i) =2*( V + pow(3, 0.5) * pot * ( 0.5));
		}
		else {
			
			ex(i) =2* pot * (3. / 2);
			}
		}
}

void Test::orthogonal(mat eigenvec) {
	//tolerance for orthogonality;
	double tol = 0.001;
	JacobiSolver jcbs;
	//dot1 = v0*v1 = 0;
	double dot1 = eigenvec(0, 0) * eigenvec(1, 0) + eigenvec(0, 1) * eigenvec(1, 1)
		+ eigenvec(0, 2) * eigenvec(1, 2);
	//dot2=v0*v0=1
	double dot2 = eigenvec(0, 0) * eigenvec(0, 0) + eigenvec(0, 1) * eigenvec(0, 1)
		+ eigenvec(0, 2) * eigenvec(0, 2);
	if (dot1 > tol) {
		cout << "Orthogonality is not preserved" << endl;
		
		
	}
	else if((dot2-1. ) > tol) {
		cout << "Orthogonality is not preserved" << endl;
		
	
	}
}


