#include "EigenvalueSolver.h"
ofstream ofile;

void EigenvalueSolver::write_to_file(string filename, vec eigen, vec ex,vec eval, int N)
{
	ofile.open(filename);
	// we write out the values to our selected file;
	ofile << "N   =   " <<   N << endl;
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "      approx:        exact:         arma:         relative error(log10)" << endl;
	for (int i = 0; i <= N - 1; i++) {
		ofile << setw(15) << setprecision(8) << eigen(i);
		ofile << setw(15) << setprecision(8) << ex(i);
		ofile << setw(15) << setprecision(8) << eval(i);
		ofile << setw(15) << setprecision(8) << log10(abs((ex(i)-eigen(i))/ex(i))) << endl;

	}
	ofile.close();
}

void EigenvalueSolver::write_to_file_2(string filename, int N, vec ex, vec wr, vec ground_values,vec armad) {
	ofile.open(filename);
	ofile << "N   =   " << N << endl;
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "             wr       approx:        exact:         arma:         relative error(log10)" << endl;
	for (int i = 0; i < 4; i++) {
		ofile << setw(15) << setprecision(8) << wr(i);
		ofile << setw(15) << setprecision(8) << ground_values(i);
		ofile << setw(15) << setprecision(8) << ex(i);
		ofile << setw(15) << setprecision(8) << armad(i);
		ofile << setw(15) << setprecision(8) << log10(abs((ex(i) - ground_values(i)) / ex(i))) << endl;
	}
	ofile.close();
}

void EigenvalueSolver::TridiagonalMatrix(mat &A, int N)
{	

	double h = 1./ (N+1);
	double hh = h * h;
	
	// For a simple code we fill in the first and last row before the loop;
	
	A(0, 0) = 2.0/hh; A(0, 1) = -1.0/hh;
	A(N - 1, N - 1) = 2.0/hh; A(N - 1, N - 2) = -1.0/hh;
	for (int i = 1; i < N - 1; i++) {
		// Here we set the diagonal elements of the matrix to the desired value;
		A(i, i) = 2.0 / hh;
		A(i, i + 1) = -1.0 / hh;
		A(i, i - 1) = -1.0 / hh;
		
	}
}
void EigenvalueSolver::QM_Matrix(mat& A, int N, vec& r, int func, double pot)
{
	double h = 50. / (N + 1);
	double hh = h * h;
	r(0) = h;

	// QM-problem for 1 electron;
	if (func == 2){
		A(0, 0) = 2.0 / hh + r(0)*r(0); A(0, 1) = -1.0 / hh;
		
		
		for (int i = 1; i < N - 1; i++) {
			// Here we set the diagonal elements of the matrix to the desired value;
			r(i) = r(i - 1) + h;;
			A(i, i) = (2.0 / hh) + (r(i)*r(i));
			A(i, i + 1) = -1.0 / hh;
			A(i, i - 1) = -1.0 / hh;
		}
		r(N - 1) = r(N - 2) + h;
		A(N - 1, N - 1) = 2.0 / hh + r(N-1)*r(N-1); A(N - 1, N - 2) = -1.0 / hh;
		
		
	}
	// QM-problem for 2 electrons;
	else if (func == 3) {
		
		A(0, 0) = (2.0 / hh) + pot * pot * r(0) * r(0) + 1. / r(0); A(0, 1) = -1.0 / hh;

		for (int i = 1; i < N - 1; i++) {
			// Here we set the diagonal elements of the matrix to the desired value;
			r(i) = r(i - 1) + h;
			A(i, i) = (2.0 / hh) + pot * pot * r(i) * r(i) + 1. / r(i);
			A(i, i + 1) = -1.0 / hh;
			A(i, i - 1) = -1.0 / hh;
		}
		r(N - 1) = r(N - 2) + h;
		A(N - 1, N - 1) = 2.0 / hh + pot* pot * r(N - 1) * r(N - 1) + 1. / r(N - 1); A(N - 1, N - 2) = -1.0 / hh;
		
	}
}