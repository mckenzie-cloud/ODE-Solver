
#include "../header/solver.h"

// A protected constructor for derived DenseOut class
Solver::Solver(std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName, size_t dim)
            : m_F{fName}, m_dim{dim}, m_yn(dim), m_X(dim), m_K1(dim), m_K2(dim), m_K3(dim), m_K4(dim), m_K5(dim), m_K6(dim), m_K7(dim), m_ynew(dim)
{}

Solver::Solver(StepSizeController::Controllers controller, std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName,
               std::vector<long double> &y0, size_t dim, long double t0, long double tFinal, long double absTol, long double relTol)
            : m_F{fName}, m_yn{y0}, m_dim{dim}, m_t{t0}, m_tFinal{tFinal}, m_absTol{absTol}, m_relTol{relTol},
              m_X(dim), m_K1(dim), m_K2(dim), m_K3(dim), m_K4(dim), m_K5(dim), m_K6(dim), m_K7(dim), m_ynew(dim), m_truncationErrors(dim), m_sci(dim)
{

    // Initialize stepsize
    m_h = initialize_stepsize(t0, y0);

    std::cout << m_h << std::endl;

    switch (controller)
    {
        // set the parameters
        case StepSizeController::STANDARD:
            set_controller_parameters(1.0L, 0.0L, 0.0L, 0.0L, 0.0L);
            break;
        case StepSizeController::H211PI:
            set_controller_parameters(1.0L / 6.0L, 1.0L / 6.0L, 0.0L, 0.0L, 0.0L);
            break;
        case StepSizeController::PI42:
            set_controller_parameters(3.0L / 5.0L, -1.0L / 5.0L, 0.0L, 0.0L, 0.0L);
            break;
        case StepSizeController::H211B: // b = 4.0L
            set_controller_parameters(1.0L / 4.0L, 1.0L / 4.0L, 0.0L, 1.0L / 4.0L, 0.0L);
            break;
        case StepSizeController::H312PID:
            set_controller_parameters(1.0L / 18.0L, 1.0L / 9.0L, 1.0L / 18.0L, 0.0L, 0.0L);
            break;
        case StepSizeController::H312B: // b = 8.0L
            set_controller_parameters(1.0L / 8.0L, 2.0L / 8.0L, 1.0L / 8.0L, 3.0L / 8.0L, 1.0L / 8.0L);
            break;
        case StepSizeController::H0321:
            set_controller_parameters(5.0L / 4.0L, 1.0L / 2.0L, -3.0L / 4.0L, -1.0L / 4.0L, -3.0L / 4.0L);
            break;
        case StepSizeController::H321:
            set_controller_parameters(1.0L / 3.0L, 1.0L / 18.0L, -5.0L / 18.0L, -5.0L / 6.0L, -1.0L / 6.0L);
            break;
    }
}

void Solver::solve()
{
    // Add the initial value
    m_tOut.push_back(m_t);
    m_yOut.push_back(m_yn);

    long double m_cerr1 = 1.0L, m_cerr2 = 1.0L, m_rh1 = 1.0L, m_rh2 = 1.0L;

    // Closed-loop system
    while (m_t < m_tFinal)
    {
        m_h = std::min(m_h, m_tFinal - m_t);

        process_doPri_steps();

        long double err_estimate_scaled = errorEstimateScaled();

        long double cerr = 1.0L / err_estimate_scaled; // inverse of the error estimates scaled

        long double rho = filter(cerr, m_cerr1, m_cerr2, m_rh1, m_rh2);

        // Save previous values for the next step.
        m_cerr2 = m_cerr1;
        m_cerr1 = cerr;

        m_rh2 = m_rh1;
        m_rh1 = rho;

        // Apply a limiter
        long double ratio = 1.0L + m_kappa * std::atan((rho - 1.0L) / m_kappa);

        if (ratio < m_acceptSF)           // Reject steps and recalculate with the new stepsize
        {
            control_stepsize(ratio);
            m_rejectedSteps++;
            continue;
        }
        else  // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.
        { 
            m_t += m_h;
            control_stepsize(ratio);
            m_acceptedSteps++;

            m_tOut.push_back(m_t);
            m_yOut.push_back(m_ynew);

            for (size_t i = 0; i < m_dim; i++)
            {
                m_yn[i] = m_ynew[i];
            }
        }
    }
}

void Solver::display_steps()
{
    std::cout << "\n\tSteps: accepted = " << m_acceptedSteps << " rejected = " << m_rejectedSteps << std::endl;
}

void Solver::display_results()
{
    std::cout << std::endl;

    constexpr auto max_precision{std::numeric_limits<long double>::digits10};

    std::cout << std::setprecision(max_precision) << std::fixed;

    for (size_t i = 0; i < m_yOut.size(); i++)
    {
        /* code */
        std::cout << "step " << i << " at t = " << m_tOut[i] << '\n';
        std::cout << "\t -> ";
        for (int j = 0; j < m_dim; j++)
        {
            std::cout << m_yOut[i][j] << " | ";
        }
        std::cout << '\n';
    }
}

void Solver::set_controller_parameters(long double b1, long double b2, long double b3, long double a2, long double a3)
{
    m_beta1 = b1 / m_k;
    m_beta2 = b2 / m_k;
    m_beta3 = b3 / m_k;
    m_alpha2 = a2;
    m_alpha3 = a3;
}

long double Solver::hairer_norm(std::vector<long double> &a, std::vector<long double> &sci)
{
    /**
     * ----- Calculate error-norm ||err|| -----
     * using the L2-Norm or the Euclidean Norm.
     * */
    long double sumOfSqrd = 0.0L;
    for (size_t i = 0; i < m_dim; i++)
    {
        sumOfSqrd = sumOfSqrd + std::pow(a[i] / sci[i], 2.0L);
    }
    return std::sqrt(sumOfSqrd / static_cast<long double>(m_dim));
}

long double Solver::initialize_stepsize(long double t0, std::vector<long double> &y0)
{
    /*
     * (a) Do one function evaluation f(t0, y0) at the initial point.
     * Then put d0 = ||y0|| and d1 = ||f(t0, y0)||, using the hairer's norm
     * with sci = Atol + |y0_i| * Rtol
     */
    std::vector<long double> f1(m_dim);
    m_F(t0, y0, f1);

    std::vector<long double> sci(m_dim);
    for (int i = 0; i < m_dim; i++)
    {
        sci[i] = m_absTol + std::fabs(y0[i]) * m_relTol;
    }

    long double d0 = hairer_norm(y0, sci);
    long double d1 = hairer_norm(f1, sci);

    // (b) As a first guess for the step size let
    long double h0 = 0.01L * (d0 / d1);

    // If either d0 or d1 is < 1e-5 we put h0 = 1e-6
    if (d0 < 1e-5 || d1 < 1e-5)
    {
        h0 = 1e-6;
    }

    // (c) Perform one explicit Euler stpe, y1 = y0 + h0 * f(t0, y0)
    std::vector<long double> y1(m_dim);
    for (int i = 0; i < m_dim; i++)
    {
        y1[i] = y0[i] + h0 * f1[i];
    }

    std::vector<long double> f2(m_dim);
    m_F(t0 + h0, y1, f2);

    /*
     * (d) Compute d2 = ||f(t0 + h0, y1) - f(t0, y0)|| / h0 as an
     * estimate of the second derivative of the solution;
     * again by using hairer's norm.
     */
    std::vector<long double> diff_f2f1(m_dim);
    for (int i = 0; i < m_dim; i++)
    {
        diff_f2f1[i] = std::fabs(f2[i] - f1[i]);
    }

    long double d2 = hairer_norm(diff_f2f1, sci) / h0;

    long double max_d1d2 = std::max(d1, d2);

    /*
     * (e) Compute a step size h1 from the relation, h1^(p+1) * max(d1, d2) = 0.0L1,
     * where p - is order of the method. If max(d1, d2) <= 10^-15,
     * put h1 = max(10^-6, h0 * 10^-3);
     */
    long double h1 = std::pow(10.0L, (-2.0L - std::log10(max_d1d2)) / m_k);
    if (max_d1d2 <= 1e-15)
    {
        h1 = std::max((long double)1e-6, h0 * 1e-3);
    }

    // f. Finally we propose as starting step size
    long double h = std::min(100.0L * h0, h1);

    return h;
}

// Dormand-Prince 5(4) steps
void Solver::process_doPri_steps()
{
    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i];
    }

    m_F(m_t, m_X, m_K1); //--------------------- 1st-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a21 * m_K1[i]);
    }

    m_F(m_t + Dopri54::c2 * m_h, m_X, m_K2); //--------------------- 2nd-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a31 * m_K1[i] + Dopri54::a32 * m_K2[i]);
    }

    m_F(m_t + Dopri54::c3 * m_h, m_X, m_K3); //--------------------- 3rd-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a41 * m_K1[i] + Dopri54::a42 * m_K2[i] +
                                  Dopri54::a43 * m_K3[i]);
    }

    m_F(m_t + Dopri54::c4 * m_h, m_X, m_K4); //--------------------- 4th-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a51 * m_K1[i] + Dopri54::a52 * m_K2[i] +
                                  Dopri54::a53 * m_K3[i] + Dopri54::a54 * m_K4[i]);
    }

    m_F(m_t + Dopri54::c5 * m_h, m_X, m_K5); //--------------------- 5th-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a61 * m_K1[i] + Dopri54::a62 * m_K2[i] +
                                  Dopri54::a63 * m_K3[i] + Dopri54::a64 * m_K4[i] + Dopri54::a65 * m_K5[i]);
    }

    m_F(m_t + m_h, m_X, m_K6); //--------------------- 6th-stage -----------------------

    // Calculate the 5th-order and 4th-order accurate solution.
    for (size_t i = 0; i < m_dim; i++)
    {
        m_ynew[i] = m_yn[i] + m_h * (Dopri54::b1 * m_K1[i] + Dopri54::b3 * m_K3[i] +
                                     Dopri54::b4 * m_K4[i] + Dopri54::b5 * m_K5[i] + Dopri54::b6 * m_K6[i]); // 5th-Order accurate solution. Used to advance the solution.
    }

    m_F(m_t + m_h, m_ynew, m_K7); //--------------------- 7th-stage -----------------------
}

long double Solver::errorEstimateScaled()
{
    // Calculate local errors
    for (size_t i = 0; i < m_dim; i++)
    {
        m_truncationErrors[i] = m_h * ((Dopri54::b1 - Dopri54::e1) * m_K1[i] + (Dopri54::b3 - Dopri54::e3) * m_K3[i] +
                                       (Dopri54::b4 - Dopri54::e4) * m_K4[i] + (Dopri54::b5 - Dopri54::e5) * m_K5[i] +
                                       (Dopri54::b6 - Dopri54::e6) * m_K6[i] -  Dopri54::e7 * m_K7[i]);
    }

    for (size_t i = 0; i < m_dim; i++)
    {
        // absTol and relTol are the desired tolerances prescribed by the user.
        m_sci[i] = m_absTol + std::max(std::fabs(m_yn[i]), std::fabs(m_ynew[i])) * m_relTol;
    }

    // local error is controlled by error-per-unit-steps (EPUS) or error-per-steps (EPS)
    // long double r = L2Norm() / m_h;  // Error-per-unit-steps (EPUS) => m_k = m_p
    long double err_estimate_scaled = hairer_norm(m_truncationErrors, m_sci); // Error-per-steps (EPS)

    return err_estimate_scaled;
}

long double Solver::filter(long double cerrPres, long double cerrOld1, long double cerrOld2, long double rho1, long double rho2)
{
    /**
     * The General controller formula for Order-Dynamics pD <= 3 with Control error filtering.
     * Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
     * https://dl.acm.org/doi/10.1145/641876.641877 -> page 22
     */
    long double result = std::pow(cerrPres, m_beta1) *
                         std::pow(cerrOld1, m_beta2) *
                         std::pow(cerrOld2, m_beta3) *
                         std::pow(rho1, -m_alpha2) *
                         std::pow(rho2, -m_alpha3);
    return result;
}

void Solver::control_stepsize(long double ratio)
{
    m_h *= ratio;
}