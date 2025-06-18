#include <iostream>

#include "../../header/solver.h"

// F(t, y)
void lorenz(long double t, std::vector<long double> &X, std::vector<long double> &Ks) {
    Ks[0] = 10.0 * (X[1] - X[0]);                // dx/dt = σ(y - x)
    Ks[1] = X[0] * (28.0 - X[2]) - X[1];         // dy/dt = x(ρ - z) - y
    Ks[2] = X[0] * X[1] - (8.0 / 3.0) * X[2];    // dz/dt = xy - βz
}

int main(void)
{
    std::vector<long double> y0 = {1.0, 1.0, 1.0};         // initial value y0
    long double t0 = 0.0;
    long double tFinal = 5.0;
    long double abstol = 1e-6, reltol = 1e-6;
    StepSizeController::Controllers controller = StepSizeController::H321; // Choose a controller
    Solver solver(controller, lorenz, y0, y0.size(), t0, tFinal, abstol, reltol);
    solver.solve();
    solver.display_results();
    solver.display_steps();

    return 0;
}