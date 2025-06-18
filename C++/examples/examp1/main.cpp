#include <iostream>

#include "../../header/solver.h"

// F(t, y)
void F(long double t, std::vector<long double> &X, std::vector<long double> &Ks)
{
    Ks[0] = X[0] - t*t + 1.0;        // y' = y - t^2 + 1
}

int main(void)
{
    std::vector<long double> y0 = {0.5};         // initial value y0
    long double t0 = 0.0;
    long double tFinal = 2.0;
    long double abstol = 1e-6, reltol = 1e-6;
    StepSizeController::Controllers controller = StepSizeController::STANDARD; // Choose a controller
    Solver solver(controller, F, y0, y0.size(), t0, tFinal, abstol, reltol);
    solver.solve();
    solver.display_results();
    solver.display_steps();

    return 0;
}