#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

const double T0 = 293;
const double g = 9.81;

double trajectory_generator(
    double t,
    double delta_t,
    double v_0,
    double theta_0,
    double D,
    double a,
    double m,
    double alpha,
    std::string filename
) {
    double x = 0, y = 0;
    double v_x = v_0 * cos(theta_0);
    double v_y = v_0 * sin(theta_0);
    double v = sqrt(pow(v_x, 2.0) + pow(v_y, 2.0));
    bool run = true;

    double x_old, y_old, t_old;

    std::ofstream trajectory(filename);

    do  {
        x_old = x;
        y_old = y;
        t_old = t;

        double F_x = -1 * D * v * v_x * pow( 1.0 - (a * y / T0), alpha);
        double F_y = -1 *  ( g * m + D * v * v_y * pow( 1.0 - (a * y / T0), alpha));

        x += delta_t*v_x;
        y += delta_t*v_y;

        v_x += delta_t*F_x/m;
        v_y += delta_t*F_y/m;
        v = sqrt(pow(v_x, 2.0) + pow(v_y, 2.0));

        t = t_old + delta_t;

        if (y < 0) {
            double r = y_old/(y_old - y);
            x = x_old + (x - x_old)*r;
            y = y_old + (y - y_old)*r;
            t = t - (1 - r)*delta_t;

            run = false;
        }
        trajectory << t << " " << x << " " << y << std::endl;
    } while (run);
    trajectory.close();

    return x;
}

int main() {
    double m = 1.0, v_0 = 100.0, theta_0 = M_PI/4.0;
    double t_max = 2 * v_0*sin(theta_0)/g;
    double x_max = v_0*cos(theta_0) * t_max;

    // trajectories with no drag force and no altitude amendement
    double n_list[] = {10.0, 20.0, 50.0, 100.0, 200.0, 500.0};

    std::ofstream file_errors("./data/no_drag_force/errors.dat");

    for (auto &n : n_list) {
        double x = trajectory_generator(0.0, t_max/n, v_0, theta_0, 0.0, 0.0, m, 0.0, 
            "./data/no_drag_force/frajectory_n" + std::to_string(n) + ".dat");
        file_errors << std::abs(x_max - x) << " " << t_max/n << std::endl;
    }

    std::cout << "entered" << std::endl;

    file_errors.close();

    // trajectories with drag force and without altitude amendement

    for (auto &D : {0.0, 1e-4, 2e-4, 5e-4, 1e-3}) {
        trajectory_generator(0.0, t_max/500.0, v_0, theta_0, D, 0.0, m, 0.0, 
            "./data/drag_force_no_altitude/frajectory_D" + std::to_string(D) + ".dat");
    }

    std::cout << "entered" << std::endl;

    // check of projectile range
    std::vector<double> angles = {};
    double angle_start = 14.0;
    while (angle_start++ < 65) angles.push_back(angle_start);

    for (auto &D : {0.0, 1e-3, 2e-3}) {
        std::ofstream file_diff_angles("./data/drag_force+angle/xmax_D" + std::to_string(D) + ".dat");
        for (auto &theta : angles) {
            double x = trajectory_generator(0.0, t_max/500.0, v_0, M_PI*theta/180.0, D, 0.0, m, 0.0, 
                "./data/drag_force+angle/unused/frajectory_D" + std::to_string(D) + "angle" + std::to_string(int(theta)) + ".dat");
            file_diff_angles << theta << " " << x << std::endl; 
        }
    }

    std::cout << "entered" << std::endl;

    // trajectories with altitude correction
    m = 20.0;
    v_0 = 700.0;
    double D = 1e-3;
    double a_list[] = {6.5*1e-3, 0.0};
    double theta_list[] {35.0, 45.0};

    for (auto &a : a_list) {
        for (auto &theta : theta_list) {
            trajectory_generator(0.0, t_max/500.0, v_0, M_PI*theta/180.0, D, a, m, 2.5,
                "./data/altitude_correction/trajectory_a" + std::to_string(a) + "theta" + std::to_string(int(theta)) + ".dat");
        }
    }

    std::cout << "entered" << std::endl;

    return 0;
}