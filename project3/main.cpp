#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#define a           0.39709         // AU
#define T_mercury   0.240846        // years
#define M_s         1.998*1e30      // kg
#define M_m         2.4*1e23        // kg
#define e           0.206           // eccentricity

#define PERIHELION_r_min        a*(1.0 - e)
#define APHELION_r_max          a*(1.0 + e)
#define PERIHELION_v_max        M_PI*2*sqrt((1.0 + e)*(1.0 + ((M_m )/ (M_s)))/(a*(1.0 - e)))
#define APHELION_v_min          M_PI*2*sqrt((1.0 - e)*(1.0 + ((M_m) / (M_s)))/(a*(1.0 + e)))


double sympletic_algorithm(
    double t_max,
    double delta_t, 
    double alpha,
    std::string dirname,
    bool plot_trajectory,
    bool plot_peri_ap_helia
) {
    // files opening
    std::ofstream trajectory;
    std::ofstream perihelia;
    std::ofstream aphelia;

    plot_trajectory ? trajectory.open(dirname+"/trajectory.dat") : trajectory.open(dirname+"unused/trajectory.dat");
    plot_peri_ap_helia ? perihelia.open(dirname+"/perihelia.dat") : perihelia.open(dirname+"unused/perihelia.dat");
    plot_peri_ap_helia ? aphelia.open(dirname+"/aphelia.dat") : aphelia.open(dirname+"unused/aphelia.dat");

    // initial conditions
    int k_max = 4;
    double x = APHELION_r_max, y = 0.0;
    double v_x = 0.0, v_y = APHELION_v_min;
    double r_0 = sqrt(x*x + y*y);
    double r_1, r_2;

    double peri_1 = 0.0, peri_2 = 0.0, time_p1 = 0.0, time_p2 = 0.0;

    double time = 0;

    // Neri parametrization
    const double a_param[] = {
        0.5 / (2.0 - pow(2.0, 1.0/3)),
        0.5 * (1.0 - pow(2.0, 1.0/3)) / (2.0 - pow(2.0, 1.0/3)),
        0.5 * (1.0 - pow(2.0, 1.0/3)) / (2.0 - pow(2.0, 1.0/3)),
        0.5 / (2.0 - pow(2.0, 1.0/3))
    };

    const double b_param[] = {
        1.0 / (2.0 - pow(2.0, 1.0/3)),
        -1.0 * pow(2.0, 1.0/3) / (2.0 - pow(2.0, 1.0/3)),
        1.0 / (2.0 - pow(2.0, 1.0/3)),
        0.0
    };

    while(time <= t_max) {
        double x_new, y_new, v_x_new, v_y_new, r_new, w_new;
        double x_old = x, y_old = y, v_x_old = v_x, v_y_old = v_y;

        for(int k = 0; k < k_max; k++) {
            x_new = x_old + a_param[k]*v_x_old*delta_t;
            y_new = y_old + a_param[k]*v_y_old*delta_t;

            r_new = sqrt(x_new*x_new + y_new*y_new);
            w_new = -4 * M_PI * M_PI * (1 + alpha / pow(r_new, 2)) / pow(r_new, 3);

            v_x_new = v_x_old + b_param[k]*w_new*x_new*delta_t;
            v_y_new = v_y_old + b_param[k]*w_new*y_new*delta_t;

            x_old = x_new;
            y_old = y_new;
            v_x_old = v_x_new;
            v_y_old = v_y_new;
        }

        if (plot_peri_ap_helia)
            if (time == 0) r_1 = sqrt(x_new*x_new + y_new*y_new);
            else {
                r_2 = sqrt(x_new*x_new + y_new*y_new);
                
                if (r_1 < r_0 && r_1 < r_2) {
                    aphelia << x << " " << y << std::endl;
                }
                else if (r_1 > r_0 && r_1 > r_2) {
                    perihelia << x << " " << y << std::endl;
                    if (peri_1 == 0.0) {
                        peri_1 = atan2(y_new, x_new);
                        time_p1 = time;
                    }
                    else if (peri_2 == 0.0) {
                        peri_2 = atan2(y_new, x_new);
                        time_p2 = time;
                    }
                }
            
                r_0 = r_1;
                r_1 = r_2;
            }

        x = x_new;
        y = y_new;
        v_x = v_x_new;
        v_y = v_y_new;

        time += delta_t;

        if(plot_trajectory) trajectory << x << " " << y << std::endl;  
    }

    trajectory.close();
    perihelia.close();
    aphelia.close();
    return (peri_2 - peri_1) / (time_p2 - time_p1); // omega(alpha)
}


int main() {
    //sympletic_algorithm(0.95*T_mercury, 1e-4, 0.0, "./data/0,95T/", true, false);
    //sympletic_algorithm(100.0*T_mercury, 1e-4, 0.0, "./data/100T/", true, false);
    //sympletic_algorithm(4.0*T_mercury, 1e-4, 0.01, "./data/alpha0,01/", true, true);

    std::ofstream prec_ang_vel("./data/prec_ang_vel/precession_angular_velocity.dat");
    double alpha_max = 0.001;
    for (int j = 0; j < 7; j++) {
        double alpha = alpha_max / pow(2.0, float(j));
        double omega = sympletic_algorithm(3, 1e-5, alpha, "./data/prec_ang_vel/", false, true);
        prec_ang_vel << alpha << " " << omega << std::endl;
    }
    prec_ang_vel.close();

    return 0;
}