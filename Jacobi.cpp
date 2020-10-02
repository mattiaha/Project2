#include "EigenvalueSolver.h"




void JacobiSolver::jsolve(int n, double conv, mat& a, vec& r, mat& v)
{
    
    cout.precision(5);
    double aip = 0, aiq = 0, vpi = 0, vqi = 0;
    double tau = 0, t = 0, s = 0, c = 0;//tan(theta), sin(theta), cos(theta)    
    int count = 1;                //count of iterations

    int p = n - 1, q = n - 2;           //off diag all same value to start
                                //pick last as first maximum
    clock_t start, end;
    double app = a(p, p);
    double aqq = a(q, q);
    double apq = a(p, q);

    start = clock();

    while (abs(apq) > conv) {
        if (count > 1) {
            apq = 0;
            find_max(a, p, q, apq, n);
        }

        //calculate sin(theta) and cos(theta)
        aqq = a(q, q);
        app = a(p, p);
        tau = (aqq - app) / (2 * apq);
        if (tau > 0)
            t = 1 / (tau + sqrt(1 + tau * tau));
        else
            t = -1 / (-tau + sqrt(1 + tau * tau));
        c = 1 / sqrt(1 + t * t);
        s = c * t;

        //calculate new matrix elements and vectors
        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                aip = a(i, p);
                aiq = a(i, q);
                a(i, p) = aip * c - aiq * s;
                a(p, i) = aip * c - aiq * s;
                a(i, q) = aiq * c + aip * s;
                a(q, i) = aiq * c + aip * s;
            }
            
            vpi = v(i, p);
            vqi = v(i, q);
            
            v(i, p) = c * vpi - s * vqi;
            v(i, q) = c * vqi + s * vpi;
        }
        a(p, p) = app * c * c - 2 * apq * c * s + aqq * s * s;
        a(q, q) = app * s * s + 2 * apq * c * s + aqq * c * c;
        a(p, q) = 0;
        a(q, p) = 0;

        count++;
    }

    end = clock();



    cout << "Diagonalization took " << count << " iterations" << endl;
    cout << scientific << "CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;

}
//get first three eigenvectors
mat JacobiSolver::get_eigenvecs(mat& a, mat& v, int n) {
    vector<double>eigenvals = get_eigenvals(a, n);
    mat vecs(3, n);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < n; j++) {
            if (a(j, j) == eigenvals[i]) {
                for (int k = 0; k < n; k++) {
                    vecs(i, k) = v(k, j);
                }
            }
        }
    }
    return vecs;
}

//get eigenvalues in order
vector<double> JacobiSolver::get_eigenvals(mat& a, int n) {
    vector<double>eigen;
    for (int i = 0; i < n; i++) {
        eigen.push_back(a(i, i));
    }
    sort(eigen.begin(), eigen.begin() + n);
    return eigen;
}



//find maximum non-diag matrix elements
void JacobiSolver::find_max(mat a, int& p, int& q, double& apq, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && abs(a(i, j)) >= abs(apq)) {
                apq = a(i, j);
                p = i;
                q = j;
            }
        }
    }
}

//print matrix and eigenvectors
void JacobiSolver::print_vals(mat A, mat v, int n, double conv) {
    cout << "A: ";
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            cout << "   ";
        }
        for (int j = 0; j < n; j++) {
            if (abs(A(i, j)) > conv)
                cout << fixed << A(i, j) << " ";
            else cout << "0.000 ";
        }
        cout << endl;
    }
    for (int i = 0; i < n; i++) {
        cout << "v" << i << ": ";
        for (int j = 0; j < n; j++) {
            if (abs(v(j, i)) > conv)
                cout << fixed << v(j, i) << " ";
            else cout << "0.000 ";
        }
        cout << endl;
    }
}

