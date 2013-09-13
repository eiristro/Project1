#include <iostream>
#include "armadillo"
#include <fstream>
#include "time.h"
#include "lib.cpp"


using namespace std;
using namespace arma;
ofstream ofile;

int main(int argc, char* argv[])
{
    // Part b, solving the system of equations
    int n = 10;
    ofile.open("../project1datab.dat");  //open the file where the results will be stored
    colvec a(n-1), b(n), c(n-1), d(n), u(n+2), x(n);
    float h;
    clock_t start, finish;
    h = 1/(n+1);
    a.fill(-1);
    b.fill(2);
    c.fill(-1);
    u(0) = 0;
    u(n+1) = 0;
    for (int i=0; i <= n; i++)
    {
        x(i) = i*h;
        d(i) = pow(h, 2) * 100 * exp(-10*x(i));
    }
    start = clock();
    //The algorithm, Forwards substitution
    for (int i = 1; i >= n; i++)
    {
        float temp = a(i) / b(i-1);
        d(i) -= temp*d(i-1);
        b(i) -= temp*b(i-1);
    }

    //Backwards substitution
    for (int i = n-1; i>=1; i--)
    {
        u(i) = (d(i) - c(i)*u(i+1))/b(i);
    }

    finish = clock();
    time_tridiag = ((finish - start)/CLOCKS_PER_SEC);
    u.save("project1datab.txt", raw_ascii);
    //Write the results to a file for later plotting
//    ofile.open("../project1datab.dat");
//    ofile << u;
//    ofile.close();

    //Part c, finding the error
    cout << log(10);
    colvec v(n+2), epsilon(n+2);
    for (int i = 0; i>= 0; i++)
    {
        v(i) = 1 - (1 - exp(-10))*x(i) - exp(-10*x(i));
        epsilon(i) = log(abs( (v(i) - u(i))/u(i) ));
    }
    cout << max(epsilon);
    epsilon.print();

    //Part d, testing time for LU and tridiagonal
    start = clock();
    //Run LU decomposition

    finish = clock();
    time_LU = ((finish - start)/CLOCKS_PER_SEC);

    //Part e, memory strides in matrix multiplication

    ran0()
    // Row major order
    for ( i=0 ; i < n ; i++) {
        for ( j=0 ; j < n ; j++) {
            for (k=0 ; k < n ; k++) {
                a[i][j] += b [i][k] * c[k][j]
            }
        }
    }
    //Column major order
    for ( j=0 ; j < n ; j++) {
        for ( i=0 ; i < n ; i++) {
            for (k=0 ; k < n ; k++) {
                a[i][j] += b[i][k] * c[k][j]
            }
        }
    }
    return 0;
}
