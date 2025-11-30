#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>


double S_temp(
    int i,
    int j,
    double delta_N,
    double S_x_min,
    double S_x_max,
    double S_y_min,
    double S_y_max,
    double S_max
    ) {
    double x = i * delta_N;
    double y = j * delta_N;

    if (x >= S_x_min && x <= S_x_max && y >= S_y_min && y <= S_y_max) {
        return S_max;
    }
    else {
        return 0.0;
    }
}

void CN_algorithm(
    double h,
    double h_omega,
    double Temp_high,
    double omega,
    std::string dirname,
    bool heater_switch
) {
    // --------------- Initial data for all runs
    const double delta_t = 10.0, T_max = 1e4;
    const double n_max = T_max / delta_t;
    const double L = 10.0;
    const int N = 51;
    const double delta_N = L / (N - 1);
    const double D = 0.1;
    const double Temp_low = 293.0, Temp_out = 273.0;
    const double S_max = 2.0;
    const int k_max = 30;
    const double TOL = 1e-8;                                                    // tolerance
    const double y_window_min = 6.0, y_window_max = 9.0;                        // window boundaries   
    const double S_x_min = 8.0, S_x_max = 8.8, S_y_min = 2.0, S_y_max = 4.0;    // heater boundaries
    const double x_sensor = 2.0, y_sensor = 8.0;                                // sensor coordinates 

    // --------------- Changable parameters
    double t = 0.0;
    double omega_old = omega;
    double **Temp = new double*[N];
    for (int i = 0; i < N; i++) {
        Temp[i] = new double[N];
        for (int j = 0; j < N; j++) {
            Temp[i][j] = Temp_out;
        }
    }

    // --------------- Output files
    std::ofstream file_Temp;
    file_Temp.open(dirname + "/tempDist.dat");

    std::cout << "Starting CN algorithm with h = " << h << " and h_omega = " << h_omega << std::endl;
    // --------------- Algorithm implementation
    for (int n = 0; n < n_max; n++) {
        std::cout << "Time step " << n + 1 << " / " << n_max << "\r";
        // calculating R_i_j
        double **R = new double*[N];
        for (int i = 0; i < N; i++) {
            R[i] = new double[N];
            for (int j = 0; j < N; j++) {  
                double T_comp = -4 * Temp[i][j];
                if (i > 0) T_comp += Temp[i - 1][j];
                if (i < N - 1) T_comp += Temp[i + 1][j];
                if (j > 0) T_comp += Temp[i][j - 1];
                if (j < N - 1) T_comp += Temp[i][j + 1];  

                R[i][j] = Temp[i][j] 
                    + D*delta_t*T_comp/(2*delta_N*delta_N) 
                    + 0.5*delta_t*omega*S_temp(i, j, delta_N, S_x_min, S_x_max, S_y_min, S_y_max, S_max);
            }
        }

        // saving old weight factor
        omega_old = omega;

        // changing the heater state (if necessary) 
        if (heater_switch == true){
            //std::cout << "Current Temp at sensor: " << Temp[int(x_sensor/delta_N)][int(y_sensor/delta_N)] << " K ";
            if (Temp[int(x_sensor/delta_N)][int(y_sensor/delta_N)] < Temp_low) {
                omega = 1;
                //std::cout << "Heater ON at time " << t << " s\n";
            }
             
            else if (Temp[int(x_sensor/delta_N)][int(y_sensor/delta_N)] > Temp_high) {
                omega = 0;
                //std::cout << "Heater OFF at time " << t << " s\n";
            }
            else {
                //std::cout << "Heater state unchanged at time " << t << " s\n";
            }
        }
        
        // Gauss-Seidel relaxation
        for (int k = 0; k < k_max; k++) {
            // updating Temp values
            for (int i = 1; i < N - 1; i++) {
                for (int j = 1; j < N - 1; j++) {
                    Temp[i][j] = ((0.5*D*delta_t/(delta_N*delta_N))*
                        (Temp[i - 1][j] + Temp[i + 1][j] + Temp[i][j - 1] + Temp[i][j + 1]) 
                        + R[i][j] + 0.5*delta_t*S_temp(i, j, delta_N, S_x_min, S_x_max, S_y_min, S_y_max, S_max))
                        / (1 + 2*D*delta_t/(delta_N*delta_N));
                }
            }
            // updating boundary conditions
            for (int i = 0; i < N - 1; i++) {
                Temp[0][i] = (h*delta_N*Temp_out + D*Temp[1][i])/(h*delta_N + D);       // west boundary
                Temp[i][N-1] = (h*delta_N*Temp_out + D*Temp[i][N-2])/(h*delta_N + D);   // north boundary
                // East boundary (x = L) - Contains the Window
                double y_curr = i * delta_N;
                double current_h = h; // Default to wall coefficient

                // Check if we are within the window coordinates
                if (y_curr >= y_window_min && y_curr <= y_window_max) {
                    current_h = h_omega; // Use window coefficient
                }
            
                // Use current_h instead of h
                Temp[N-1][i] = (current_h*delta_N*Temp_out + D*Temp[N-2][i])/(current_h*delta_N + D);  
                Temp[i][0] = (h*delta_N*Temp_out + D*Temp[i][1])/(h*delta_N + D);       // south boundary
            }

            // upgrading value at corners
            Temp[0][0] = 0.5*(Temp[0][1] + Temp[1][0]);
            Temp[0][N-1] = 0.5*(Temp[0][N-2] + Temp[1][N-1]);
            Temp[N-1][0] = 0.5*(Temp[N-2][0] + Temp[N-1][1]);
            Temp[N-1][N-1] = 0.5*(Temp[N-2][N-1] + Temp[N-1][N-2]); 

            // calculating the residual
            double c = 0.0;
            for (int i = 1; i < N - 2; i++) {
                for (int j = 1; j < N - 2; j++) {
                    double T_comp = -4 * Temp[i][j];
                    if (i > 0) T_comp += Temp[i - 1][j];
                    if (i < N - 1) T_comp += Temp[i + 1][j];
                    if (j > 0) T_comp += Temp[i][j - 1];
                    if (j < N - 1) T_comp += Temp[i][j + 1]; 

                    double L = Temp[i][j] 
                        - 0.5*D*delta_t*T_comp/(delta_N*delta_N)
                        - 0.5*delta_t*omega*S_temp(i, j, delta_N, S_x_min, S_x_max, S_y_min, S_y_max, S_max);
                    c += std::pow(L - R[i][j], 2)*delta_N*delta_N;
                }
            }

            // checking the convergence
            if (std::sqrt(c) < TOL) break;
        }

        // calculating supplied and lost (through window) energy
        double E_supplied = 0.0;
        double E_lost = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                E_supplied += omega*S_temp(i, j, delta_N, S_x_min, S_x_max, S_y_min, S_y_max, S_max);
            }
            if (i*delta_N >= y_window_min && i*delta_N <= y_window_max) {
                E_lost += (Temp[N-1][i] - Temp_out);
            }
        }
        E_supplied *= 0.5*(omega + omega_old)*(delta_N*delta_N)*delta_t;
        E_lost *= delta_N*delta_t*h_omega;

        // changing time
        t += delta_t;

        // saving temp data
        file_Temp << t << " "  << Temp[int(x_sensor/delta_N)][int(y_sensor/delta_N)] << " " << E_supplied << " " << E_lost << "\n";
        
        switch (int(t)) {
            case 2500:
            case 3000:
            case 3500:
            case 4000:
                {
                    std::ofstream file_tempDist;
                    file_tempDist.open(dirname + "/tempDist_" + std::to_string(int(t)) + ".dat");
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            file_tempDist << i*delta_N << " " << j*delta_N << " " << Temp[i][j] << std::endl;
                        }
                    }
                    file_tempDist.close();
                }
                break;
        }
        // freeing the space
        for (int i = 0; i < N; i++) 
            delete[] R[i];
    }

    // closing files
    file_Temp.close();

    // freeing the space
    for (int i = 0; i < N; i++)
        delete[] Temp[i];
}

int main() {
    // // task 1
    // CN_algorithm(0.0, 0.0, 1e4, 1, "./data/1", true);

    // // task 2
    // CN_algorithm(0.002, 0.002, 1e4, 1, "./data/2", true);

    // // task 3
    // CN_algorithm(0.002, 0.03, 1e4, 1, "./data/3", true);
    
    // task 4
    CN_algorithm(0.002, 0.03, 298, 1, "./data/4", true);

    // task 5
    CN_algorithm(0.0, 1.0, 298, 1, "./data/5", true);
    
    // task 6
    CN_algorithm(0.0, 1.0, 1e4, 1, "./data/6", true);
    return 0;
}