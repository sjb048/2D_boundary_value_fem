//g++ -I /opt/homebrew/opt/eigen/include/eigen3 possion.cpp -o possion
//./poisson

#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;

// Function to create the matrix A
MatrixXd createMatrixA(int n, double h) {
    MatrixXd A = MatrixXd::Zero(n + 1, n + 1);
    for (int i = 0; i <= n; ++i) {
        if (i > 0) A(i, i - 1) = -1 / h;   // Sub-diagonal
        A(i, i) = 2 / h;                   // Diagonal
        if (i < n) A(i, i + 1) = -1 / h;   // Super-diagonal
    }
    return A;
}

// Function to create the matrix R
MatrixXd createMatrixR(int n, double kappa) {
    MatrixXd R = MatrixXd::Zero(n + 1, n + 1);
    R(0, 0) = kappa;      // R(0,0) for left boundary
    R(n, n) = kappa;      // R(n,n) for right boundary
    return R;
}

// Function to create the vector b
VectorXd createVectorB(int n, double h) {
    VectorXd b = VectorXd::Zero(n + 1);
    b(0) = h / 2;         // b(0)
    b(n) = h / 2;         // b(n)
    for (int i = 1; i < n; ++i) {
        b(i) = h;         // All other entries are h (since f = 1)
    }
    return b;
}

// Function to create the vector r
VectorXd createVectorR(int n, double g0, double gL, double kappa) {
    VectorXd r = VectorXd::Zero(n + 1);
    r(0) = kappa * g0;    // Contribution from left boundary
    r(n) = kappa * gL;    // Contribution from right boundary
    return r;
}

// Function to solve the boundary value problem
void solvePoisson(double L, int n, double g0, double gL, double kappa) {
    double h = L / n;  // Element size

    // Create matrices and vectors
    MatrixXd A = createMatrixA(n, h);
    MatrixXd R = createMatrixR(n, kappa);
    VectorXd b = createVectorB(n, h);
    VectorXd r = createVectorR(n, g0, gL, kappa);

    // Solve (A + R)ζ = b + r
    MatrixXd M = A + R;                          // M = A + R
    VectorXd zeta = M.colPivHouseholderQr().solve(b + r); // ζ = (A + R)⁻¹(b + r)

    // Output coefficients ζ
    cout << "Coefficients ζ for n = " << n << ":" << endl;
    for (int i = 0; i <= n; ++i) {
        cout << "ζ[" << i << "] = " << zeta(i) << endl;
    }

    // Save coefficients to a CSV file for plotting
    ofstream outFile("coefficients_n" + to_string(n) + ".csv");
    outFile << "Index,Zeta\n";
    for (int i = 0; i <= n; ++i) {
        outFile << i << "," << zeta(i) << "\n";
    }
    outFile.close();
}

int main() {
    double L = 1.0;          // Domain length
    double g0 = 0.25;       // Boundary condition at x = 0
    double gL = 0.0;        // Boundary condition at x = L
    double kappa = 1e6;     // Constant for Dirichlet conditions

    // Test for different values of n
    for (int n = 5; n <= 50; n += 5) {
        solvePoisson(L, n, g0, gL, kappa);
    }

    return 0;
}

