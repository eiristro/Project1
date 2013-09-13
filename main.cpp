#include <iostream>
#include <iomanip>
#include "armadillo"
#include <fstream>
#include "time.h"
#include "lib.cpp"
#include <math.h>
#include <typeinfo>


using namespace std;
using namespace arma;
ofstream ofile;

int main()
{
    // Part b, solving the system of equations
    int n = 10;
    ofile.open("../epsilon2.txt");
    while (n<=1e5)
    {
        // Declaring variables for the tridiagonal algorithm
        // and time measurement, the right hand side is called d
        colvec a(n+1), b(n+1), c(n+1), d(n+2), v(n+2), x(n+2);
        float h;
        clock_t start, finish;

        // Initializing values
        h = 1.0/(n+1);
        a.fill(-1);
        b.fill(2);
        c.fill(-1);
        v(0) = 0;
        v(n+1) = 0;

        for (int i=0; i < n+2; i++)  // Calculating x and d
        {
            x(i) = i*h;
            d(i) = pow(h, 2) * 100 * exp(-10*x(i));
        }
        start = clock();

        //The algorithm, Forwards substitution,  starting at i = 2
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
        finish = clock();

        // Calculating the time the algorithm ran
        float time_tridiag = ((finish - start)/(float)CLOCKS_PER_SEC);

//        //Write the results to a file for later plotting
//        char buffer [50];
//        sprintf(buffer, "../project1datab%d.txt", n);
//        v.save(buffer, raw_ascii);

        //Part c, finding the error
        colvec u(n+2), epsilon(n+2);
        for (int i = 0; i < n+2; i++)
        {
            u(i) = 1 - (1 - exp(-10))*x(i) - exp(-10*x(i));
            epsilon(i) = log10(abs( (v(i) - u(i))/u(i) ));
        }

        // Output the greatest relative error, v(0) and v(n+1) is by definition 0
        // so the relative error for these terms is always zero
        float max_epsilon = max(epsilon.rows(2, n));
        cout << max_epsilon << endl;
//        ofile << max_epsilon << endl;

        //Part d, testing time for LU and tridiagonal
        //Setting up the matrix, the tridiagonal algorithm only used the vectors
        double ** A;
        A = new double * [n];
        for (int i = 0; i < n; i++) {
            A[i] = new double[n];
        }

        // Assigning values to the matrix: 2 on the diagonal, -1 just off it
        // and 0 everywhere else
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

        // Declaring the variables ludcmp and lubksb need
        int indx[n];
        double dd[n];
        // bb has the same values as d, but lubsb does not accept armadillo vectors
        double bb[n];
        for (int i = 0; i<n; i++) {
            bb[i] = pow(h, 2) * 100 * exp(-10*(i+1)*h);
        }
        start = clock();

        // Run LU decomposition
        ludcmp(A, n, indx, dd);
        lubksb(A, n, indx, bb);

        finish = clock();

        // Calculating the time the algorithm ran
        float time_LU = ((finish - start)/(float)CLOCKS_PER_SEC);

        // Freeing memory
        for (int i = 0; i < n; i++) {
            delete[] A[i];
        }
        delete[] A;

        // Updating n and outputting the times
        n = n*10;
        cout << time_tridiag << " Tridiagonal, ";
        cout << time_LU << " LU, " << endl;
    }

    //Part e, memory strides in matrix multiplication
    int n = 1e3;
    clock_t start, finish;

    // Declaring the matrices
    mat B(n,n), C(n,n), D(n,n);

    // Initializing with random values
    long int q = 132;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B(i, j) = ran0(&q);
            C(i, j) = ran0(&q);
        }
    }

    // Row major order
    start = clock();
    for (int i=0 ; i < n ; i++) {
        for (int j=0 ; j < n ; j++) {
            for (int k=0 ; k < n ; k++) {
                D(i,j) += B(i, k) * C(k, j);
            }
        }
    }
    finish = clock();
    float time_RowMajor = ((float(finish - start))/CLOCKS_PER_SEC);

    //Column major order
    start = clock();
    for (int j=0 ; j < n ; j++) {
        for (int i=0 ; i < n ; i++) {
            for (int k=0 ; k < n ; k++) {
                D(i,j) += B(i, k) * C(k, j);
            }
        }
    }

    finish = clock();
    float time_ColumnMajor = ((finish - start)/(float)CLOCKS_PER_SEC);

    cout << time_RowMajor << " Row major, ";
    cout << time_ColumnMajor << " Column major" << endl;

    ofile.close();
    return 0;
}
