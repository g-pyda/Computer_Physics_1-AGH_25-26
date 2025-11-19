#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

void save_potential(
    std::string dirname,
    double **v,
    double **rho,
    int N,
    double delta,
    double L
) {
    std::ofstream potential_file(dirname + "/V(x,y).dat");
    std::ofstream error_file(dirname + "/err(x,y).dat");
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            potential_file << x * delta - L << " " << y * delta - L << " " << v[y][x] << std::endl;
            error_file << x * delta - L << " " << y * delta - L << " " << (delta * delta * v[y][x]) + rho[y][x] << " " << (delta * delta * v[y][x]) << " " << rho[y][x]  << std::endl;
        }
    }
    potential_file.close();
    error_file.close();
}

void save_SOR_iterations(
    std::string dirname,
    int iter,
    double s,
    double delta_s
) {
    std::ofstream SOR_iterations_file(dirname + "/SOR_iterations.dat", std::ios_base::app);
    SOR_iterations_file << iter << " " << s << " " << delta_s << std::endl;
    SOR_iterations_file.close();
}

double stopping_criterion(
    double **v,
    double **rho,
    int N,
    double delta
) {
    double s = 0.0;
    for (int y = 1; y < N-1; y++) {
        for (int x = 1; x < N-1; x++) {
            double Ex = -1 * (v[y][x+1] - v[y][x-1]) / (2 * delta);
            double Ey = -1 * (v[y+1][x] - v[y-1][x]) / (2 * delta);
            s += ((0.5 * (Ex*Ex + Ey*Ey)) - v[y][x] * rho[y][x]) * delta * delta;
        }
    }
    return s;
}

void SOR(
    double omega,
    double k1,
    double k2,
    double k3,
    double k4,
    double A_rho,
    bool Neumann_boundary,
    bool Dirichlet_boundary,
    std::string dirname

) {
    // initial data
    const int L = 4, 
        N = 100;
    double eps = 1.0, 
        TOL = 1e-8, 
        Kmax = 10000, 
        s_new = 0.0, 
        s_old = 1.0,
        delta = (double)L * 2.0 / (N - 1);
    double **v = new double*[N], 
        **rho = new double*[N];
    for (int y = 0; y < N; y++) {
        v[y] = new double[N];
        rho[y] = new double[N];
        for (int x = 0; x < N; x++) {
            v[y][x] = 0.0;
            rho[y][x] = A_rho * (x*delta - L) * (y*delta - L) * exp(-1 * (pow(x*delta - L, 2) + pow(y*delta - L, 2)));
        }
    }

    // dirichlet boundary conditions - if applicable
    if (Dirichlet_boundary == true) {
        for (int i = 0; i < N; i++) {
            v[0][i] = sin(k4 * (i * delta) * M_PI / (2 * L));
            v[N-1][i] = sin(k2 * (i * delta) * M_PI / (2 * L));
            v[i][0] = sin(k1 * (i * delta) * M_PI / (2 * L));
            v[i][N-1] = sin(k3 * (i * delta) * M_PI / (2 * L));
        }
    }

    // relaxation loop
    for (int k = 1; k <= Kmax; k++) {
        std::cout << "SOR iteration: " << k << std::endl;
        // Neumann boundary conditions - if applicable
        if (Neumann_boundary == true)
            for (int i = 1; i < N-1; i++) 
                v[i][N-1] = v[i][N-2];
        
        // upgrade of the solution
        for (int y = 1; y < N-1; y++)
            for (int x = 1; x < N-1; x++) 
                v[y][x] = (1 - omega) * v[y][x] + (omega / 4.0) * (v[y+1][x] + v[y-1][x] + v[y][x+1] + v[y][x-1] + rho[y][x] * delta * delta / eps);

        // check of covergence
        s_new = stopping_criterion(v, rho, N, delta);
        double delta_s = fabs(s_new - s_old) / fabs(s_old);
        save_SOR_iterations(dirname, k, s_new, delta_s);
        if (delta_s < TOL && k > 1) 
            break;
        s_old = s_new;
    }

    save_potential(dirname, (double **)v, (double **)rho, N, delta, L);
}

int main() {
    SOR(1.5, 1.0, -1.0, 1.0, -1.0, 0.0, false, true, "./data/initial_run");
    SOR(1.5, 1.0, -1.0, 1.0, -1.0, 0.0, true, true, "./data/initial_run_with_boundaries");
    // double omegas[] = {1.0, 1.3, 1.6, 1.9};
    // for (double omega : omegas) {
    //     std::string dirname = "./data/omega_" + std::to_string(omega);
    //     SOR(omega, 1.0, -1.0, 1.0, -1.0, 0.0, true, true, dirname);
    // }
    // SOR(1.5, 0.0, 0.0, 0.0, 0.0, 1.0, false, true, "./data/homogeneous_bc");
    // SOR(1.5, 1.0, -1.0, 1.0, -1.0, 1.0, false, true, "./data/nohomogeneous_bc");
    return 0;
}