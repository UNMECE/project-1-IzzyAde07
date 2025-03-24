#include <iostream>
#include <cmath>

// structure for the capacitor values
struct _capacitor {
    double *time;       // time array
    double *voltage;    // voltage array
    double *current;    // current array
    double C;           // capacitance value
};
typedef struct _capacitor Capacitor;

// Function for constant current charging
void simulateConstantCurrent(Capacitor &cap, double I, double dt, int numSteps) {
    // Initializing the arrays
    cap.time = new double[numSteps];
    cap.voltage = new double[numSteps];
    cap.current = new double[numSteps];

    // Initial conditions
    cap.time[0] = 0.0;
    cap.voltage[0] = 0.0; // Voltage at t=0 is 0
    cap.current[0] = I;   // Constant current

    // for loop to iterate over number of time steps
    for (int i = 1; i < numSteps; i++) {
        cap.time[i] = cap.time[i - 1] + dt;
        cap.voltage[i] = cap.voltage[i - 1] + cap.current[i - 1] * dt / cap.C;
        cap.current[i] = I; // Current remains constant
    }
}

// Function for constant voltage charging
void simulateConstantVoltage(Capacitor &cap, double V0, double R, double dt, int numSteps) {
    // Initializing the arrays
    cap.time = new double[numSteps];
    cap.voltage = new double[numSteps];
    cap.current = new double[numSteps];

    // Initial conditions
    cap.time[0] = 0.0;
    cap.voltage[0] = 0.0; // Voltage at t=0 is 0
    cap.current[0] = V0 / R; // Initial current

    // for loop to iterate over number of time steps
    for (int i = 1; i < numSteps; i++) {
        cap.time[i] = cap.time[i - 1] + dt;
        cap.current[i] = cap.current[i - 1] - (cap.current[i - 1] / (R * cap.C)) * dt;
        cap.voltage[i] = cap.voltage[i - 1] + cap.current[i - 1] * dt / cap.C;
    }
}

// Function prints results every 200 timesteps
void printResults(const Capacitor &cap, int numSteps) {
    for (int i = 0; i < numSteps; i += 200) {
        std::cout << "Time: " << cap.time[i] << " s, "
                  << "Voltage: " << cap.voltage[i] << " V, "
                  << "Current: " << cap.current[i] << " A" << std::endl;
    }
}

int main() {
    // variables for the RCI equation
    double dt = 1e-10;          // Time step
    double finalTime = 5e-6;    // Final time
    int numSteps = finalTime / dt; // Number of timesteps
    double R = 1e3;             // Resistance 
    double C = 100e-12;         // Capacitance 
    double I = 1e-2;            // Constant current 
    double V0 = 10.0;           // Constant voltage 

    // Capacitor struct
    Capacitor cap;
    cap.C = C;

    // constant current charging
    std::cout << "Constant Current Charging:" << std::endl;
    simulateConstantCurrent(cap, I, dt, numSteps);
    printResults(cap, numSteps);

    // Free memory
    delete[] cap.time;
    delete[] cap.voltage;
    delete[] cap.current;

    // constant voltage charging
    std::cout << "\nConstant Voltage Charging:" << std::endl;
    simulateConstantVoltage(cap, V0, R, dt, numSteps);
    printResults(cap, numSteps);

    // Free memory
    delete[] cap.time;
    delete[] cap.voltage;
    delete[] cap.current;

    return 0;
}