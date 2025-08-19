#include <iostream>
#include <algorithm>

#include "../../header/solver.h"
#include "../../header/dense_out.h"

// F(t, y)
void F(long double t, std::vector<long double> &X, std::vector<long double> &Ks)
{
    Ks[0] = X[0] - t*t + 1.0L;        // y' = y - t^2 + 1
}

int main(void)
{
    std::vector<long double> y0 = {0.5L};         // initial value y0
    long double t0 = 0.0L;
    long double tFinal = 2.0L;
    long double abstol = 1e-6, reltol = 1e-6;
    StepSizeController::Controllers controller = StepSizeController::STANDARD; // Choose a controller
    Solver solver(controller, F, y0, y0.size(), t0, tFinal, abstol, reltol);
    solver.solve();
    solver.display_results();
    solver.display_steps();
    
    /*
     * Using the dense out class to "query" the solution at any time (t) in the integration interval.
     * if t_query = t_i => yn(t_i). Otherwise, if t_query = t_{i+1} => yn(t_{t+1}). Where yn(t_i) is the solution at time (t_i).
    */
    // long double t_query = 1.7115918395774204035L;
    // long double t_query = 2.0L;
    // DenseOut ds(F, solver.m_yOut, solver.m_tOut, y0.size());
    // std::vector<long double> solP = ds.dense_eval_sol(t_query);
    // ds.printDebug();

    // std::cout << "Solution: " << "\n";

    // for (size_t i = 0; i < solP.size(); i++)
    // {
    //     std::cout << solP[i] << '\n';
    // }

    return 0;
}