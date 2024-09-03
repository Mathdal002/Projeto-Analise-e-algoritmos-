#include <iostream>
#include <vector>

using namespace std;

// Função para decompor a matriz A em L e U
bool luDecomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        // Parte de U
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }
        // Parte de L
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal principal de L
            else {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
    return true;
}

// Função para resolver L*y = b
vector<double> forwardSubstitution(vector<vector<double>>& L, vector<double>& b) {
    int n = b.size();
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
    return y;
}

// Função para resolver U*x = y
vector<double> backwardSubstitution(vector<vector<double>>& U, vector<double>& y) {
    int n = y.size();
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}

int main() {
    vector<vector<double>> A = {
        {4, -1, 0, 0, 0, 0, 0, 0, 0, 0},
        {-1, 4, -1, 0, 0, 0, 0, 0, 0, 0},
        {0, -1, 4, -1, 0, 0, 0, 0, 0, 0},
        {-1, 0, 0, 4, -1, 0, 0, 0, 0, 0},
        {0, -1, 0, -1, 4, -1, 0, 0, 0, 0},
        {0, 0, -1, 0, -1, 4, -1, 0, 0, 0},
        {0, 0, 0, -1, 0, -1, 4, -1, 0, 0},
        {0, 0, 0, 0, -1, 0, -1, 4, -1, 0},
        {0, 0, 0, 0, 0, -1, 0, -1, 4, -1},
        {0, 0, 0, 0, 0, 0, -1, 0, -1, 4}
    };
    vector<double> b = {-110, -30, -40, -110, 0, -15, -90, -25, -55, -65};

    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    if (!luDecomposition(A, L, U)) {
        cout << "Decomposição LU falhou." << endl;
        return 1;
    }

    vector<double> y = forwardSubstitution(L, b);
    vector<double> x = backwardSubstitution(U, y);

    cout << "Solução usando Fatoração LU:" << endl;
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << endl;

    return 0;
}