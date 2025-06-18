
#include "../header/solver.h"

Solver::Solver(StepSizeController::Controllers controller, std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName,
               std::vector<long double> &y0, size_t dim, long double t0, long double tFinal, long double absTol, long double relTol, bool denseOut)
            : m_F{fName}, m_yn{y0}, m_dim{dim}, m_t{t0}, m_tFinal{tFinal}, m_absTol{absTol}, m_relTol{relTol}, m_denseOut{denseOut},
              m_X(dim), m_K1(dim), m_K2(dim), m_K3(dim), m_K4(dim), m_K5(dim), m_K6(dim), m_K7(dim), m_ynew(dim), m_truncationErrors(dim), m_sci(dim)
{

    // Initialize stepsize
    m_h = initialize_stepsize(t0, y0);

    std::cout << m_h << std::endl;

    switch (controller)
    {
        // set the parameters
        case StepSizeController::STANDARD:
            set_controller_parameters(1.0, 0.0, 0.0, 0.0, 0.0);
            break;
        case StepSizeController::H211PI:
            set_controller_parameters(1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0);
            break;
        case StepSizeController::PI42:
            set_controller_parameters(3.0 / 5.0, -1.0 / 5.0, 0.0, 0.0, 0.0);
            break;
        case StepSizeController::H211B: // b = 4.0
            set_controller_parameters(1.0 / 4.0, 1.0 / 4.0, 0.0, 1.0 / 4.0, 0.0);
            break;
        case StepSizeController::H312PID:
            set_controller_parameters(1.0 / 18.0, 1.0 / 9.0, 1.0 / 18.0, 0.0, 0.0);
            break;
        case StepSizeController::H312B: // b = 8.0
            set_controller_parameters(1.0 / 8.0, 2.0 / 8.0, 1.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0);
            break;
        case StepSizeController::H0321:
            set_controller_parameters(5.0 / 4.0, 1.0 / 2.0, -3.0 / 4.0, -1.0 / 4.0, -3.0 / 4.0);
            break;
        case StepSizeController::H321:
            set_controller_parameters(1.0 / 3.0, 1.0 / 18.0, -5.0 / 18.0, -5.0 / 6.0, -1.0 / 6.0);
            break;
    }
}

void Solver::solve()
{
    // Add the initial value
    m_tOut.push_back(m_t);
    m_yOut.push_back(m_yn);

    long double m_cerr1 = 1.0, m_cerr2 = 1.0, m_rh1 = 1.0, m_rh2 = 1.0;

    // Closed-loop system
    while (m_t < m_tFinal)
    {
        m_h = std::min(m_h, m_tFinal - m_t);

        long double err_estimate_scaled = process();

        long double cerr = 1.0 / err_estimate_scaled; // inverse of the error estimates scaled

        long double rho = filter(cerr, m_cerr1, m_cerr2, m_rh1, m_rh2);

        // Save previous values for the next step.
        m_cerr2 = m_cerr1;
        m_cerr1 = cerr;

        m_rh2 = m_rh1;
        m_rh1 = rho;

        // Apply a limiter
        long double ratio = 1.0 + m_kappa * std::atan((rho - 1.0) / m_kappa);

        if (ratio < m_acceptSF)
        {
            // Reject steps and recalculate with the new stepsize
            control_stepsize(ratio);
            m_rejectedSteps++;
            continue;
        }
        else
        { 
            // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.
            if (m_denseOut) // do 1-more extra steps at the middle between yn and yn1.
            {
                /*
                 * theta ∈ [0, 1], theta = 0 => yn, theta = 1 => yn1
                 * interpolate at open-interval theta ∈ (0, 1)
                 * un+1(t + theta*h) = yn + h * sum(bi(theta)*Ki), i = 1...s, theta ∈ (0, 1)
                 */
                long double theta = 0.5; 
                std::vector<long double> extraSteps = interpolate(theta, m_h);
                m_tOut.push_back(m_t + theta * m_h);
                m_yOut.push_back(extraSteps);
            }

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

    constexpr auto max_precision{std::numeric_limits<long double>::digits10 + 1};

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
    long double sumOfSqrd = 0.0;
    for (size_t i = 0; i < m_dim; i++)
    {
        sumOfSqrd = sumOfSqrd + std::pow(a[i] / sci[i], 2.0);
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
    long double h0 = 0.01 * (d0 / d1);

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
     * (e) Compute a step size h1 from the relation, h1^(p+1) * max(d1, d2) = 0.01,
     * where p - is order of the method. If max(d1, d2) <= 10^-15,
     * put h1 = max(10^-6, h0 * 10^-3);
     */
    long double h1 = std::pow(10.0, (-2.0 - std::log10(max_d1d2)) / m_k);
    if (max_d1d2 <= 1e-15)
    {
        h1 = std::max((long double)1e-6, h0 * 1e-3);
    }

    // f. Finally we propose as starting step size
    long double h = std::min(100.0 * h0, h1);

    return h;
}

long double Solver::process()
{

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i];
    }

    m_F(m_t, m_X, m_K1); //--------------------- 1ST-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a21 * m_K1[i]);
    }

    m_F(m_t + Dopri54::c2 * m_h, m_X, m_K2); //--------------------- 2ND-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a31 * m_K1[i] + Dopri54::a32 * m_K2[i]);
    }

    m_F(m_t + Dopri54::c3 * m_h, m_X, m_K3); //--------------------- 3RD-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a41 * m_K1[i] + Dopri54::a42 * m_K2[i] +
                                  Dopri54::a43 * m_K3[i]);
    }

    m_F(m_t + Dopri54::c4 * m_h, m_X, m_K4); //--------------------- 4TH-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a51 * m_K1[i] + Dopri54::a52 * m_K2[i] +
                                  Dopri54::a53 * m_K3[i] + Dopri54::a54 * m_K4[i]);
    }

    m_F(m_t + Dopri54::c5 * m_h, m_X, m_K5); //--------------------- 5TH-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        m_X[i] = m_yn[i] + m_h * (Dopri54::a61 * m_K1[i] + Dopri54::a62 * m_K2[i] +
                                  Dopri54::a63 * m_K3[i] + Dopri54::a64 * m_K4[i] + Dopri54::a65 * m_K5[i]);
    }

    m_F(m_t + m_h, m_X, m_K6); //--------------------- 6TH-stage -----------------------

    // Calculate the 5th-order and 4th-order accurate solution.
    for (size_t i = 0; i < m_dim; i++)
    {
        m_ynew[i] = m_yn[i] + m_h * (Dopri54::b1 * m_K1[i] + Dopri54::b3 * m_K3[i] +
                                     Dopri54::b4 * m_K4[i] + Dopri54::b5 * m_K5[i] + Dopri54::b6 * m_K6[i]); // 5th-Order accurate solution. Used to advance the solution.
    }

    m_F(m_t + m_h, m_ynew, m_K7); //--------------------- 7TH-stage -----------------------

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

std::vector<long double> Solver::interpolate(long double theta, long double hPresent)
{

    const long double C1 = 5.0 * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
    const long double C3 = 100.0 * (882725551.0 - 15701508.0 * theta) / 32700410799.0;
    const long double C4 = 25.0 * (443332067.0 - 31403016.0 * theta) / 1880347072.0;
    const long double C5 = 32805.0 * (23143187.0 - 3489224.0 * theta) / 199316789632.0;
    const long double C6 = 55.0 * (29972135.0 - 7076736.0 * theta) / 822651844.0;
    const long double C7 = 10.0 * (7414447.0 - 829305.0 * theta) / 29380423.0;

    long double theta_sqr = std::pow(theta, 2.0);
    long double term1 = theta_sqr * (3.0 - 2.0 * theta);
    long double term2 = theta_sqr * std::pow(theta - 1.0, 2.0);
    long double term3 = theta * std::pow(theta - 1.0, 2.0);
    long double term4 = (theta - 1.0) * std::pow(theta, 2.0);

    long double b1Theta = term1 * Dopri54::b1 + term3 - term2 * C1;
    long double b3Theta = term1 * Dopri54::b3 + term2 * C3;
    long double b4Theta = term1 * Dopri54::b4 - term2 * C4;
    long double b5Theta = term1 * Dopri54::b5 + term2 * C5;
    long double b6Theta = term1 * Dopri54::b6 - term2 * C6;
    long double b7Theta = term4 + term2 * C7;

    std::vector<long double> solution(m_dim);
    for (size_t i = 0; i < m_dim; i++)
    {
        /* code */
        solution[i] = m_yn[i] + hPresent * (b1Theta * m_K1[i] + b3Theta * m_K3[i] + b4Theta * m_K4[i] +
                                            b5Theta * m_K5[i] + b6Theta * m_K6[i] + b7Theta * m_K7[i]);
    }
    return solution;
}