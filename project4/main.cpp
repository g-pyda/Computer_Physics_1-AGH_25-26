#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

std::vector<double> trapezoidal(
    double delta_t,
    std::vector<double> N,
    std::vector<double> lambda,
    double tol
) {
    std::vector<double> N_new = N, delta_N(3);
    double eps_max;
    int it_max = 50, it = -1;

    const double G[3][3] = {
        {1.0 + 0.5*delta_t*lambda[0], 0.0, 0.0},
        {-0.5*delta_t*lambda[0], 1.0 + 0.5*delta_t*lambda[1], 0.0},
        {0.0, -0.5*delta_t*lambda[1], 1.0 + 0.5*delta_t*lambda[2]}
    };

    const double det_G = G[0][0]*(G[1][1]*G[2][2] - G[1][2]*G[2][1]) - G[0][1]*(G[1][0]*G[2][2] - G[1][2]*G[2][0]) + G[0][2]*(G[1][0]*G[2][1] - G[1][1]*G[2][0]);

    do {
        it++;
        std::vector<double> F(3);
        F[0] = N_new[0] - N[0] + 0.5 * delta_t * lambda[0] * (N[0] + N_new[0]);
        F[1] = N_new[1] - N[1] - 0.5 * delta_t * ((lambda[0]*N[0] - lambda[1]*N[1]) + (lambda[0]*N_new[0] - lambda[1]*N_new[1]));
        F[2] = N_new[2] - N[2] - 0.5 * delta_t * ((lambda[1]*N[1] - lambda[2]*N[2]) + (lambda[1]*N_new[1] - lambda[2]*N_new[2]));

        delta_N[0] = -1.0 * F[0] * (G[1][1]*G[2][2] - G[1][2]*G[2][1]) / det_G
                     + F[1] * (G[0][1]*G[2][2] - G[0][2]*G[2][1]) / det_G
                     - F[2] * (G[0][1]*G[1][2] - G[0][2]*G[1][1]) / det_G;
        delta_N[1] = F[0] * (G[1][0]*G[2][2] - G[1][2]*G[2][0]) / det_G
                     -1.0 * F[1] * (G[0][0]*G[2][2] - G[0][2]*G[2][0]) / det_G
                     + F[2] * (G[0][0]*G[1][2] - G[0][2]*G[1][0]) / det_G;
        delta_N[2] = -1.0 * F[0] * (G[1][0]*G[2][1] - G[1][1]*G[2][0]) / det_G
                     + F[1] * (G[0][0]*G[2][1] - G[0][1]*G[2][0]) / det_G
                     - F[2] * (G[0][0]*G[1][1] - G[0][1]*G[1][0]) / det_G;  

        for(int i = 0; i < 3; i++) N_new[i] += delta_N[i];

        if (fabs(delta_N[0]) >= fabs(delta_N[1]) && fabs(delta_N[0]) >= fabs(delta_N[2]))
            eps_max = fabs(delta_N[0]);
        else if (fabs(delta_N[1]) >= fabs(delta_N[0]) && fabs(delta_N[1]) >= fabs(delta_N[2]))
            eps_max = fabs(delta_N[1]);
        else
            eps_max = fabs(delta_N[2]);
    } while (eps_max > tol && it < it_max);

    for (int i = 0; i < 3; i++) N[i] = N_new[i];

    return N;
}

void adaptive_timestep(
    double tolerance_trapezoidal,
    double tolerance_timestep,
    double t_max,
    double delta_t,
    std::vector<double> N,
    std::vector<double> lambda,
    std::string N_directory,
    std::string timestep_directory
) {
    double S = 0.9, TOL = tolerance_timestep, t = 0.0;

    std::vector<double> eps(3), N_1(3), N_2(3), N_2_2(3);
    double eps_max, delta_t_old;

    std::ofstream N_file(N_directory);
    std::ofstream timestep_file(timestep_directory);

    while ( t < t_max) {
        do {
            N_1 = trapezoidal(delta_t, N, lambda, tolerance_trapezoidal);
            N_2_2 = trapezoidal(delta_t / 2.0, N, lambda, tolerance_trapezoidal);
            N_2 = trapezoidal(delta_t / 2.0, N_2_2, lambda, tolerance_trapezoidal);
            
            for(int i = 0; i < 3; i++) eps[i] = (N_2[i] - N_1[i]) / 3.0;

            if (fabs(eps[0]) >= fabs(eps[1]) && fabs(eps[0]) >= fabs(eps[2]))
                eps_max = fabs(eps[0]);
            else if (fabs(eps[1]) >= fabs(eps[0]) && fabs(eps[1]) >= fabs(eps[2]))
                eps_max = fabs(eps[1]);
            else
                eps_max = fabs(eps[2]);

            delta_t_old = delta_t;
            delta_t *= S * pow(TOL / eps_max, 1.0 / 3.0);

        } while (TOL / eps_max < 1);

        for(int i = 0; i < 3; i++) N[i] = N_2[i] + eps[i];
        
        t += delta_t_old;

        N_file << t << " " << N[0] << " " << N[1] << " " << N[2] << std::endl;
        timestep_file << t << " " << delta_t_old << std::endl;
    }
    N_file.close();
    timestep_file.close();
}

void gen_exact_solution(
    std::string N_directory,
    double t_max,
    double delta_t,
    std::vector<double> N,
    std::vector<double> lambda
) {
    std::ofstream N_exact(N_directory);
    for (double t = 0.0; t <= t_max; t += 5*delta_t) {
        double N0_exact = N[0] * exp(-lambda[0]*t);

        double N1_exact = -1.0 * N[0] * lambda[0] * exp(-lambda[0]*t) 
                            / (lambda[0] - lambda[1]) 
                         + (lambda[0] * (N[0] + N[1]) - lambda[1]*N[1]) * exp(-lambda[1]*t) 
                            / (lambda[0] - lambda[1]);

        double N2_exact = -1.0 * N[0] * lambda[0] * lambda[1] * exp(-lambda[0]*t) 
                            / ((lambda[0] - lambda[1]) * (lambda[2] - lambda[0]))  
                        + lambda[1] * (lambda[0]*(N[0] + N[1]) - lambda[1]*N[1]) * exp(-lambda[1]*t) 
                            / ((lambda[0] - lambda[1]) * (lambda[2] - lambda[1]))
                        + (lambda[0]*lambda[1]*(N[0] + N[1] + N[2]) - lambda[0]*lambda[2]*N[2] - lambda[1]*lambda[2]*(N[1] + N[2]) + pow(lambda[2],2.0)*N[2]) * exp(-lambda[2]*t)
                            / ((lambda[2] - lambda[1]) * (lambda[2] - lambda[0]));

        N_exact << t << " " << N0_exact << " " << N1_exact << " " << N2_exact << std::endl;
    }
    N_exact.close();
}

int main() {
    // initial conditions
    std::vector<double> N = {1.0, 0.0, 0.0}, lambda(3);
    double delta_t = 1e-2, t_max = 200.0, tol_trap = 1e-10, tol_timestep;

    //test of correctness
    lambda = {1.0, 5.0, 50.0};
    tol_timestep = 1e-4;

    adaptive_timestep(tol_trap, tol_timestep, t_max, delta_t, N, lambda, "data/0/N.dat", "data/0/timestep.dat");
    gen_exact_solution("data/0/N_exact.dat", t_max, delta_t, N, lambda);

    //repetition 1
    lambda = {100.0, 1.0, 0.01};
    tol_timestep = 1e-6;

    adaptive_timestep(tol_trap, tol_timestep, t_max, delta_t, N, lambda, "data/1/N.dat", "data/1/timestep.dat");
    gen_exact_solution("data/1/N_exact.dat", t_max, delta_t, N, lambda);

    //repetition 2
    tol_timestep = 1e-3;

    adaptive_timestep(tol_trap, tol_timestep, t_max, delta_t, N, lambda, "data/2/N.dat", "data/2/timestep.dat");
    gen_exact_solution("data/2/N_exact.dat", t_max, delta_t, N, lambda);

    return 0;
}