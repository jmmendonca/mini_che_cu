#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cstring>

struct MiniChemParams {
    double T_in;
    double P_in;
    double t_step;
    int n_step;
    int n_sp;
    char data_file[200];
    char sp_file[200];
    char network[200];
    char net_dir[200];
    char met[200];
};

struct MiniChemVMRParams {
    int CE_IC;
    char IC_file[200];
    double* VMR_IC;
};

void miniChDlsode(double T_in, double P_in, double t_step, double* VMR, char* network) {
    // Implementation of the mini_ch_dlsode function
    // ... (you need to implement this function)
}










int main() {
    std::ifstream file("mini_chem.nml");
    MiniChemParams params;
    MiniChemVMRParams VMRParams;

    //params.network = "OH";
    //params.network = "CHO";
    std::strcpy(params.network, "NCHO");

    params.T_in = 1500.0;
    params.P_in = 1.0e6;

    params.t_step = 60;
    params.n_step = 500000;

    params.n_sp = 12;

    std::strcpy(params.data_file, "chem_data/mini_chem_data_NCHO.txt");
    std::strcpy(params.sp_file, "chem_data/mini_chem_sp_NCHO.txt");
    std::strcpy(params.net_dir, "chem_data/1x/");
    std::strcpy(params.met, "1x");

    VMRParams.CE_IC = 1;
    std::strcpy(VMRParams.IC_file, "chem_data/IC/mini_chem_IC_FastChem_1x.txt");

    //VMRParams.VMR_IC = {0.0, 0.8, 0.0, 0.0, 0.2};
    //VMRParams.VMR_IC = {0.0, 0.8, 0.0, 0.0, 0.1, 0.0, 0.1, 0.0, 0.0};
    
    // Allocate memory for the array
    VMRParams.VMR_IC = new double[params.n_sp];

    // Initialize the array with values
    double initialVMR[] = {0.0, 0.9975, 0.001074, 0.0, 0.0, 0.0, 0.0, 0.00059024, 0.0, 0.00014159, 0.0, 0.0};

    // Copy the values into the array
    std::copy(std::begin(initialVMR), std::end(initialVMR), VMRParams.VMR_IC);

    std::cout << "T [K], P [bar], t_step, n_step, n_sp :" << std::endl;
    std::cout << params.T_in << " " << params.P_in << " " << params.t_step << " " << params.n_step << " " << params.n_sp << std::endl;

    // Initial time
    double t_now = 0.0;

    // Read the reaction and species list
    // (you need to implement the read_react_list function)

    // Save the initial conditions to file
    // Rescale IC to 1
    double sum_VMR_IC = std::accumulate(VMRParams.VMR_IC, VMRParams.VMR_IC + params.n_sp, 0.0);
    for (int i = 0; i < params.n_sp; ++i)
        VMRParams.VMR_IC[i] /= sum_VMR_IC;

    std::cout << "IC: ";
    for (int i = 0; i < params.n_sp; ++i)
        std::cout << VMRParams.VMR_IC[i] << " ";

    std::cout << sum_VMR_IC << std::endl;

    // Give initial conditions to VMR array
    double* VMR = new double[params.n_sp];
    std::memcpy(VMR, VMRParams.VMR_IC, params.n_sp * sizeof(double));

    // Do time marching loop
    // - this loop emulates what a call to the model is like in the GCM
    for (int n = 1; n <= params.n_step; ++n) {
        // Update time
        t_now += params.t_step;

        // Time now
        std::cout << n << " " << params.n_step << " " << t_now << std::endl;

        // Scale VMR to 1
        double sum_VMR = std::accumulate(VMR, VMR + params.n_sp, 0.0);
        for (int i = 0; i < params.n_sp; ++i)
            VMR[i] /= sum_VMR;

        // Call dlsode - bdf method
        miniChDlsode(params.T_in, params.P_in, params.t_step, VMR, params.network);
        std::cout << "dlsode: ";
        for (int i = 0; i < params.n_sp; ++i)
            std::cout << VMR[i] << " ";

        std::cout << sum_VMR << std::endl;

        // Scale VMR to 1
        for (int i = 0; i < params.n_sp; ++i)
            VMR[i] /= sum_VMR;
    }

    // Clean up
    delete[] VMR;
    delete[] VMRParams.VMR_IC;

    return 0;
}
