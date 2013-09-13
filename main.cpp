#include <iostream>
#include <iomanip>
#include "armadillo"
#include <fstream>
#include "time.h"
#include "lib.cpp"
#include <math.h>


using namespace std;
using namespace arma;
ofstream ofile;

int main()
{
    // Part b, solving the system of equations
    int n = 10;
    ofile.open("../epsilon2.txt");
//    while (n<=1e2) {
    ofile.open("../project1datab.dat");  //open the file where the results will be stored
    colvec a(n+1), b(n+1), c(n+1), d(n+2), v(n+2), x(n+2);
    float h;
    clock_t start, finish;
    h = 1.0/(n+1);
    a.fill(-1);
    b.fill(2);
    c.fill(-1);
    v(0) = 0;
    v(n+1) = 0;
    for (int i=0; i < n+2; i++)
    {
        x(i) = i*h;
        d(i) = pow(h, 2) * 100 * exp(-10*x(i));
    }
    clock_t t;
    t = clock();
    //The algorithm, Forwards substitution
    for (int i = 2; i < n+1; i++)
    {
        float temp = a(i) / b(i-1);
        d(i) -= temp*d(i-1);
        b(i) -= temp*c(i-1);
    }
    //Backwards substitution
    for (int i = n; i>0; i--)
    {
        v(i) = (d(i) - c(i)*v(i+1))/b(i);
    }
    t = clock() - t;
    cout << t << ", " << ((float)t)/CLOCKS_PER_SEC;
    //Write the results to a file for later plotting
    char buffer [50];
    sprintf(buffer, "../project1datab%d.txt", n);
    v.save(buffer, raw_ascii);

    //Part c, finding the error
    colvec u(n+2), epsilon(n+2);
    for (int i = 0; i < n+2; i++)
    {
        u(i) = 1 - (1 - exp(-10))*x(i) - exp(-10*x(i));
        epsilon(i) = log10(abs( (v(i) - u(i))/u(i) ));
    }

//    cout << max(epsilon) << endl;
//    ofile << max(epsilon);

    n = n*10;

    ofile.close();

    //Part d, testing time for LU and tridiagonal
    //Setting up the matrix, I only used the vectors earlier
    double ** A;
    A = new double * [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
    }
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            if (j == i) {
                A[i][j] = 2;
            }
            else if ((j == i+1) || (j == i-1)) {
                A[i][j] = -1;
            }
            else {
                A[i][j] = 0;
            }
        }
    }
    int indx[n];
    double dd[n];
    double *bb;
    bb = new double [n];
    for (int i = 0; i<n; i++) {
        bb[i] = pow(h, 2) * 100 * exp(-10*(i+1)*h);
    }
    start = clock();
    //Run LU decomposition
    ludcmp(A, n, indx, dd);
//    lubksb(A, n, indx, bb);

    finish = clock();
    float time_LU = ((finish - start)/CLOCKS_PER_SEC);
    cout << time_LU << endl;
    for (int i = 0; i < n; i++) {
        delete[] A[i];
    }
    delete[] A;
//    }

    //Part e, memory strides in matrix multiplication

//    mat B(4,4), C(4,4), D(4,4); //Using D because A is already in use
//    long int* q,* z;
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; i < n; i++) {
//            B(i, j) = ran0(q);
//            C(i, j) = ran0(z);
//        }
//    }
//    // Row major order
//    start = clock();
//    for (int i=0 ; i < n ; i++) {
//        for (int j=0 ; j < n ; j++) {
//            for (int k=0 ; k < n ; k++) {
//                D(i,j) += B(i, k) * C(k, j);
//            }
//        }
//    }
//    finish = clock();
//    float time_RowMajor = ((finish - start)/CLOCKS_PER_SEC);
//    //Column major order
//    start = clock();
//    for (int j=0 ; j < n ; j++) {
//        for (int i=0 ; i < n ; i++) {
//            for (int k=0 ; k < n ; k++) {
//                D(i,j) += B(i, k) * C(k, j);
//            }
//        }
//    }
//    finish = clock();
//    float time_ColumnMajor = ((finish - start)/CLOCKS_PER_SEC);
//    ((finish - start)/CLOCKS_PER_SEC);
    ofile.close();
    return 0;
}
