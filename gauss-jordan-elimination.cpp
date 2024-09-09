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
    int n = a.size();
    int m = b[0].size();
    vector<int> indxc(n), indxr(n), ipiv(n);

    for (int j = 0; j < n; j++) ipiv[j] = 0;

    for (int i = 0; i < n; i++) {
        cout << "Starting iteration " << i + 1 << " of " << n << endl;
        double big = 0.0;
        int irow = 0, icol = 0;

        for (int j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }

        cout << "Pivot element found at (" << irow << ", " << icol << ") with value " << big << endl;

        ++(ipiv[icol]);

        if (irow != icol) {
            for (int l = 0; l < n; l++) std::swap(a[irow][l], a[icol][l]);
            for (int l = 0; l < m; l++) std::swap(b[irow][l], b[icol][l]);
        }

        indxr[i] = irow;
        indxc[i] = icol;

        if (a[icol][icol] == 0.0) {
            cerr << "Error: Singular matrix encountered" << endl;
            throw runtime_error("gaussj: Singular matrix");
        }

        double pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;

        for (int l = 0; l < n; l++) a[icol][l] *= pivinv;
        for (int l = 0; l < m; l++) b[icol][l] *= pivinv;

        for (int ll = 0; ll < n; ll++) {
            if (ll != icol) {
                double dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (int l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
                for (int l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
            }
        }

        cout << "Iteration " << i + 1 << " completed" << endl;
    }

    for (int l = n - 1; l >= 0; l--) {
        if (indxr[l] != indxc[l]) {
            for (int k = 0; k < n; k++) {
                std::swap(a[k][indxr[l]], a[k][indxc[l]]);
            }
        }
    }

    cout << "Gauss-Jordan elimination completed successfully." << endl;
}

int main() {
    try {
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

        cout << "\nStarting Gauss-Jordan elimination..." << endl;
        gaussj(A, B);

        // Print transformed matrices
        cout << "\nTransformed A (should be identity matrix):" << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << A[i][j] << " ";
            }
            cout << endl;
        }

        cout << "\nTransformed B (solution):" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i+1 << " = " << B[i][0] << endl;
        }

        cout << "\nProgram completed successfully." << endl;
    }
    catch (const exception& e) {
        cerr << "An error occurred: " << e.what() << endl;
    }

    cout << "Press Enter to exit...";
    cin.get();

    return 0;
}



