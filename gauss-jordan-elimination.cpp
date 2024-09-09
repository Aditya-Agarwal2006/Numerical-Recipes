//2.1 in book

#include <cmath>
#include "nrx/NR.H"
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;

void gaussj(vector<vector<double>>& a, vector<vector<double>>& b)
// this is for solving linear equations by gauss-jordan elimination. 
{
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv;

    int n=a.size();
    int m=b[0].size(); // Changed from b.size() to b[0].size()
    vector<int> indxc(n),indxr(n),ipiv(n); //int arrays used for pivotign
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {  //main loop over columns to be reduced
        big=0.0;
        for(j=0;j<n;j++) { //outer loop of the search for pivot element
            if (ipiv[j] != 1) 
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {  // Changed from a[j][k] to a(j,k)
                            big=fabs(a[j][k]);       // Changed from a[j][k] to a(j,k)
                            irow=j;
                            icol=k;
                        }
                    }
                }
        }
    ++(ipiv[icol]);
    
    if (irow != icol) {
        for (l=0;l<n;l++) {
            std::swap(a[irow][l], a[icol][l]); // Changed from SWAP to std::swap
            std::swap(b[irow][l], b[icol][l]); // Changed from SWAP to std::swap
        }
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) throw runtime_error("gaussj: singular matrix");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<n;ll++) 
        if (ll != icol) {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
        }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                std::swap(a[k][indxr[l]], a[k][indxc[l]]); // Changed from SWAP to std::swap
    }
}

int main() {
    const int n = 3; // Size of the square matrix

    // Create test matrices
    vector<vector<double>> A(n, vector<double>(n));
    vector<vector<double>> B(n, vector<double>(1));

    // Initialize A (coefficient matrix)
    A[0] = {3, 2, -1};
    A[1] = {2, -2, 4};
    A[2] = {-1, 0.5, -1};

    // Initialize B (constant terms)
    B[0] = {1};
    B[1] = {-2};
    B[2] = {0};

    // Print original matrices
    cout << "Original A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nOriginal B:" << endl;
    for (int i = 0; i < n; i++) {
        cout << B[i][0] << endl;
    }

    // Solve the system
    gaussj(A, B);

    // Print results
    cout << "\nSolution:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i+1 << " = " << B[i][0] << endl;
    }

    return 0;
}



