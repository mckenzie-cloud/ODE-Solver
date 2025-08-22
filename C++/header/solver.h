
#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <limits>

#include "../header/dopri54.h"
#include "../header/stepsizeController.h"


class Solver
{

protected:
    std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> m_F; // f(t, y)
    std::vector<long double> m_yn, m_X, m_K1, m_K2, m_K3, m_K4, m_K5, m_K6, m_K7, m_ynew;
    long double m_h, m_t;
    size_t m_dim;
    Solver(std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName, size_t dim);  // constructor for derived DenseOut class
    void process_doPri_steps();

private:
    const long double m_p = 4.0;         // the order corresponding to the RK method
    const long double m_k = m_p + 1.0;   // EPS => k = p + 1 and EPUS => k = p
    const long double m_kappa = 1.0;     // kappa ∈ [0.7, 2] as suggested in the literature
    const long double m_acceptSF = 0.81; // accept safety factor

    std::vector<long double> m_truncationErrors, m_sci;

    long double m_tFinal, m_absTol, m_relTol;

    long double m_beta1, m_beta2, m_beta3, m_alpha2, m_alpha3;

    int m_acceptedSteps = 0;
    int m_rejectedSteps = 0;

private:
    void set_controller_parameters(long double b1, long double b2, long double b3, long double a2, long double a3);

    long double initialize_stepsize(long double t0, std::vector<long double> &y0);

    long double hairer_norm(std::vector<long double> &a, std::vector<long double> &sci);

    long double errorEstimateScaled();

    long double filter(long double cerrPres, long double cerrOld1, long double cerrOld2, long double rho1, long double rho2);

    void control_stepsize(long double ratio);


public:
    std::vector<std::vector<long double>> m_yOut;            // accumulate solution steps
    std::vector<long double> m_tOut;                         // accumulate time steps

    // Solver(controllerType, f(t, y), y0, t0, tFinal, absTolerance=1E-6, relTolerance=1E-4, denseOut=false)
    Solver(StepSizeController::Controllers controller, std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName,
           std::vector<long double> &y0, size_t dim, long double t0, long double tFinal, long double absTol = 1e-6, long double relTol = 1e-4);

    void solve();
    void display_steps();
    void display_results();
};

#endif