//2.1 in book

#include <cmath>
#include "nrx/NR.H"
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;

void gaussj(vector<vector<double>>& a, vector<vector<double>>& b) {
    int n = a.size();
    int m = b[0].size();

    for (int i = 0; i < n; i++) {
        // Find pivot
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(a[k][i]) > abs(a[max_row][i])) {
                max_row = k;
            }
        }

        // Swap maximum row with current row
        for (int k = i; k < n; k++) {
            swap(a[max_row][k], a[i][k]);
        }
        for (int k = 0; k < m; k++) {
            swap(b[max_row][k], b[i][k]);
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double c = -a[k][i] / a[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    a[k][j] = 0;
                } else {
                    a[k][j] += c * a[i][j];
                }
            }
            for (int j = 0; j < m; j++) {
                b[k][j] += c * b[i][j];
            }
        }
    }

    // Solve equation Ax=b using back substitution
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < m; j++) {
            b[i][j] /= a[i][i];
        }
        for (int k = i - 1; k >= 0; k--) {
            for (int j = 0; j < m; j++) {
                b[k][j] -= a[k][i] * b[i][j];
            }
        }
    }

    // Make the diagonal elements 1
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                a[i][j] = 0;
            } else {
                a[i][j] = 1;
            }
        }
    }
}

int main() {
    try {
        int n;
        cout << "Enter the size of the square matrix: ";
        cin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<vector<double>> B(n, vector<double>(1));

        cout << "Enter the elements of matrix A (" << n << "x" << n << "):" << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << "A[" << i << "][" << j << "]: ";
                cin >> A[i][j];
            }
        }

        cout << "Enter the elements of vector B (" << n << " elements):" << endl;
        for (int i = 0; i < n; i++) {
            cout << "B[" << i << "]: ";
            cin >> B[i][0];
        }

        cout << "\nOriginal A:" << endl;
        for (const auto& row : A) {
            for (double val : row) cout << val << " ";
            cout << endl;
        }

        cout << "\nOriginal B:" << endl;
        for (const auto& row : B) cout << row[0] << endl;

        gaussj(A, B);

        cout << "\nTransformed A (should be identity matrix):" << endl;
        for (const auto& row : A) {
            for (double val : row) cout << val << " ";
            cout << endl;
        }

        cout << "\nTransformed B (solution):" << endl;
        for (int i = 0; i < n; i++) cout << "x" << i+1 << " = " << B[i][0] << endl;

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    cout << "Press Enter to exit...";
    cin.get();
    cin.get();  // Add an extra cin.get() to catch the Enter key press

    return 0;
}



