#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

std::tuple<double, double> f_ODE(
    std::tuple<double,double>u, 
    double k, 
    double m,
    double alpha,
    double F_0,
    double Omega_ext,
    double t
){
    double v = std::get<1>(u);
    double F_tot_m = -1 * ( k * std::get<0>(u) + alpha * std::get<1>(u)) / m + F_0 * sin(Omega_ext*t);

    return std::tuple<double,double>{v, F_tot_m};
}

double runge_kutta_4(
    double alpha,
    double F_ext,
    double Omega_ext,
    double t_max,
    double N,
    std::string directory,
    bool energy,
    bool x_t,
    bool v_t,
    bool x_v,
    bool amplitude
){
    // initial data
    double delta_t = t_max / N;
    double k = 1.0, m = 1.0, omega_0 = 1.0, x_0 = 1.0, v_0 = 0.0, t = 0.0;
    double x_1, x_2, x_max;
    std::tuple<double, double> u{x_0, v_0};

    double time;
    std::tuple<double,double> omega, k_1, k_2, k_3, k_4;

    // opening the save files
    std::ofstream energy_file(directory+"energy.dat");
    std::ofstream v_t_file(directory+"v(t).dat");
    std::ofstream x_t_file(directory+"x(t).dat");
    std::ofstream x_v_file(directory+"x(v).dat");
    std::ofstream amplitude_file(directory+"ampl.dat");

    if (energy) {
        double E_kin = m*pow(v_0, 2)/2,
            E_pot = k*pow(x_0, 2)/2;
        energy_file << t << " " << E_kin << " " << E_pot << " " << E_kin + E_pot << std::endl;
    }
    if (v_t) v_t_file << t << " " << v_0 << std::endl;
    if (x_t) x_t_file << t << " " << x_0 << std::endl;
    if (x_v) x_v_file << v_0 << " " << x_0 << std::endl;

    // algorithm
    for(int i = 0; i < N; i++) {
        // k_1
        omega = u;
        time = t;
        k_1 = f_ODE(omega, k, m, alpha, F_ext, Omega_ext, time);
        //std::cout << std::get<0>(k_1) << " " << std::get<1>(k_1) <<std::endl;

        // k_2
        omega = std::tuple<double,double>{
            std::get<0>(u) + delta_t*std::get<0>(k_1) / 2,
            std::get<1>(u) + delta_t*std::get<1>(k_1) / 2
        };
        time = t + delta_t/2;
        k_2 = f_ODE(omega, k, m, alpha, F_ext, Omega_ext, time);
        //std::cout << std::get<0>(k_2) << " " << std::get<1>(k_2) <<std::endl;

        // k_3
        omega = std::tuple<double,double>{
            std::get<0>(u) + delta_t*std::get<0>(k_2) / 2,
            std::get<1>(u) + delta_t*std::get<1>(k_2) / 2
        };
        time = t + delta_t/2;
        k_3 = f_ODE(omega, k, m, alpha, F_ext, Omega_ext, time);
        //std::cout << std::get<0>(k_3) << " " << std::get<1>(k_3) <<std::endl;

        // k_4
        omega = std::tuple<double,double>{
            std::get<0>(u) + delta_t*std::get<0>(k_3),
            std::get<1>(u) + delta_t*std::get<1>(k_3)
        };
        time = t + delta_t;
        k_4 = f_ODE(omega, k, m, alpha, F_ext, Omega_ext, time);
        //std::cout << std::get<0>(k_4) << " " << std::get<1>(k_4) <<std::endl;

        // update of solution and time
        u = std::tuple<double,double>{
            std::get<0>(u) + delta_t*(std::get<0>(k_1) + 2*std::get<0>(k_2) + 2*std::get<0>(k_3) + std::get<0>(k_4))/6,
            std::get<1>(u) + delta_t*(std::get<1>(k_1) + 2*std::get<1>(k_2) + 2*std::get<1>(k_3) + std::get<1>(k_4))/6
        };
        t += delta_t;
        //std::cout << std::get<0>(u) << " " << std::get<1>(u) <<std::endl;

        // printing the data to files
        if (energy) energy_file << t << " " << m*pow(std::get<1>(u), 2)/2 << " " << k*pow(std::get<0>(u), 2)/2 << " " << m*pow(std::get<1>(u), 2)/2 + k*pow(std::get<0>(u), 2)/2 << std::endl;
        if (v_t) v_t_file << t << " " << std::get<1>(u) << std::endl;
        if (x_t) x_t_file << t << " " << std::get<0>(u) << std::endl;
        if (x_v) x_v_file << std::get<1>(u) << " " << std::get<0>(u) << std::endl; 

        // calculating overall maximum
        if (i == 0) x_1 = std::get<0>(u);
        else {
            x_2 = std::get<0>(u);
            if (x_1 > x_0 && x_1 > x_2) x_max = x_1;
            x_0 = x_1;
            x_1 = x_2;
        }
    }

    // closing files
    if (energy) energy_file.close();
    if (v_t) v_t_file.close();
    if (x_t) x_t_file.close();
    if (x_v) x_v_file.close();

    //std::cout << "Fuction done" << std::endl;

    return x_max;
}

int main() {
    double x_0 = 1.0, v_0 = 0.0, t_max = 50, N = 1e2, alpha = 0.1, m = 1.0, omega_0 = 1.0;
    double delta_t = t_max / N, tau =2*m/alpha;
    double omega_d = sqrt(omega_0*omega_0 - 1 / pow(tau,2));
    double phase = atan((x_0 + v_0*tau) / (x_0 * omega_d * tau));

    // // TESTING CORRECTNESS
    // runge_kutta_4(0.0, 0.0, 0.0, 50, 1e4, "./data/testing_correctness/", true, true, true, true, false);
    // // -- different alpha
    // runge_kutta_4(0.1, 0.0, 0.0, 50, 1e4, "./data/testing_correctness-alpha0.1/", true, true, false, true, false);
    // // ---- generating exact solution
    std::ofstream exact_file("./data/testing_correctness-alpha0.1/exact.dat");
    for (double t = 0.0; t <= t_max; t += delta_t) exact_file << t << " "
           << sqrt(pow(x_0,2) + pow((x_0 + v_0*tau)/(omega_d*tau),2))
              * exp(-t/tau)
              * cos(omega_d*t + phase)
           << std::endl;
    // // ENERGY DISSIPATION
    // runge_kutta_4(1e-4, 0.0, 0.0, 50, 1e4, "./data/energy_dissipation-alpha/0.0001/", false, true, false, false, false);
    // runge_kutta_4(0.1, 0.0, 0.0, 50, 1e4, "./data/energy_dissipation-alpha/0.1/", false, true, false, false, false);
    // runge_kutta_4(0.5, 0.0, 0.0, 50, 1e4, "./data/energy_dissipation-alpha/0.5/", false, true, false, false, false);
    // runge_kutta_4(1.95, 0.0, 0.0, 50, 1e4, "./data/energy_dissipation-alpha/1.95/", false, true, false, false, false);

    // // EXTERNAL PERIODIC FORCE
    // runge_kutta_4(1.0, 1.0, 0.5, 1e3, 2e5, "./data/external_force/1.0/", false, true, false, false, false);

    // std::ofstream xmax_alpha_0_01("./data/external_force/alpha0.01.dat");
    // std::ofstream xmax_alpha_0_1("./data/external_force/alpha0.1.dat");
    // std::ofstream xmax_alpha_0_5("./data/external_force/alpha0.5.dat");
    // std::ofstream xmax_alpha_1_0("./data/external_force/alpha1.0.dat");

    // for(double omega = 0.1; omega <= 2.0; omega += 0.01) {
    //     xmax_alpha_1_0 << omega << " " << runge_kutta_4(1.0, 1.0, omega, 1e3, 2e5, "./data/external_force-alpha/1.0/", false, false, false, false, false) << std::endl;
    //     xmax_alpha_0_01 << omega << " " << runge_kutta_4(0.01, 1.0, omega, 1e3, 2e5, "./data/external_force-alpha/0.01/", false, false, false, false, false) << std::endl;
    //     xmax_alpha_0_1 << omega << " " << runge_kutta_4(0.1, 1.0, omega, 1e3, 2e5, "./data/external_force-alpha/0.1/", false, false, false, false, false) << std::endl;
    //     xmax_alpha_0_5 << omega << " " << runge_kutta_4(0.5, 1.0, omega, 1e3, 2e5, "./data/external_force-alpha/0.5/", false, false, false, false, false) << std::endl;
    //     std::cout << omega << " done" << std::endl;
    // }

    // xmax_alpha_0_01.close();
    // xmax_alpha_0_1.close();
    // xmax_alpha_0_5.close();
    // xmax_alpha_1_0.close();

    return 0;
}