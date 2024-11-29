#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

void FTCS(double alpha, double u, double dx, double dt, int nx, int nt) {
    vector<double> T_old(nx, 0.0), T_new(nx, 0.0);

    // Initial condition: Linear variation from 0 to 100 along x
    for (int i = 0; i < nx; i++) {
        T_old[i] = (100.0 / (nx - 1)) * i;
    }

    // Boundary conditions
    T_old[0] = 0.0;
    T_old[nx - 1] = 100.0;

    // Open file to write results
    ofstream file("diffusion_output.txt");

    // Time stepping
    for (int t = 0; t < nt; t++) {
        // Apply FTCS scheme
        for (int i = 1; i < nx - 1; i++) {
            T_new[i] = T_old[i] - (u * dt / (2 * dx)) * (T_old[i + 1] - T_old[i - 1]) +
                       (alpha * dt / (dx * dx)) * (T_old[i + 1] - 2 * T_old[i] + T_old[i - 1]);
        }

        // Enforce boundary conditions
        T_new[0] = 0.0;
        T_new[nx - 1] = 100.0;

        // Update for the next time step
        T_old = T_new;

        // Write results at specific time steps
        if (t == 10 || t == 20 || t == 100) { // Adjust based on required time steps
            file << "Time step: " << t << endl;
            for (int i = 0; i < nx; i++) {
                file << "x = " << i * dx << ", T = " << T_new[i] << endl;
            }
            file << endl;
        }
    }

    file.close();
}

int main() {
    double alpha = 0.01;
    double u = 0.05;
    double dx = 0.1;
    double dt = 0.5;
    int nx = static_cast<int>(1.0 / dx) + 1;
    int nt = static_cast<int>(50.0 / dt);

    FTCS(alpha, u, dx, dt, nx, nt);

    cout << "Results written to 'diffusion_output.txt'" << endl;
    return 0;
}
