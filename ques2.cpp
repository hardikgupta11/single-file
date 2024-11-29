#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

void solveSecondOrder(double B, double dx, int nx) {
    vector<double> T(nx, 0.0), analytical(nx, 0.0);
    double A1 = 1.0, A2 = 1.0;

    // Matrix coefficients
    double r = B * B * dx * dx;

    // Boundary conditions
    T[0] = 1.0;
    T[nx - 1] = 1.0;

    // Numerical solution using finite difference method
    for (int i = 1; i < nx - 1; i++) {
        T[i] = 0.5 * (T[i - 1] + T[i + 1]) / (1 + r); // Iterative FD solver
    }

    // Analytical solution
    for (int i = 0; i < nx; i++) {
        double x = i * dx;
        analytical[i] = A1 * cos(B * x) + A2 * sin(B * x);
    }

    // Write results to file
    ofstream file("second_order_output.txt");
    file << "x\tNumerical\tAnalytical" << endl;
    for (int i = 0; i < nx; i++) {
        file << i * dx << "\t" << T[i] << "\t" << analytical[i] << endl;
    }
    file.close();
}

int main() {
    double B = M_PI / 2.0;
    double dx = 0.1;
    int nx = static_cast<int>(1.0 / dx) + 1;

    solveSecondOrder(B, dx, nx);

    cout << "Results written to 'second_order_output.txt'" << endl;
    return 0;
}
