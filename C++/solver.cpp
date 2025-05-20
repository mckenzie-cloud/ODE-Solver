#include <iostream>
#include <cstdlib>
#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>

/*
 *
 * Copyright (c) 2024 Mckenzie Regalado - All Rights Reserved.
 * Philippines

 * References :
 * [1] J.R. Dormand, & P.J. Prince (1980). A family of embedded Runge-Kutta formulae. Journal of Computational and Applied Mathematics, 6(1), 19-26.
 * [2] Gustafsson, K. (1991). Control Theoretic Techniques for Stepsize Selection in Explicit Runge-Kutta Methods. ACM Transactions on Mathematical Software, 17(4), 533-554.
 * [3] Söderlind, G. (2003). Digital filters in adaptive time-stepping. ACM Trans. Math. Softw., 29(1), 1–26.
 * [4] Söderlind, G., & Wang, L. (2006). Adaptive Time-Stepping and Computational Stability. Journal of Computational and Applied Mathematics, 185, 225-243. https://doi.org/10.1016/j.cam.2005.03.008
 * [5] Hairer, Ernst et. al. 𝘚𝘰𝘭𝘷𝘪𝘯𝘨 𝘖𝘳𝘥𝘪𝘯𝘢𝘳𝘺 𝘋𝘪𝘧𝘧𝘦𝘳𝘦𝘯𝘵𝘪𝘢𝘭 𝘌𝘲𝘶𝘢𝘵𝘪𝘰𝘯𝘴 𝘐: 𝘕𝘰𝘯𝘴𝘵𝘪𝘧𝘧 𝘗𝘳𝘰𝘣𝘭𝘦𝘮𝘴, Springer Science & Business Media, 1987,
 * "Automatic Step Size Control" page-167, "Dense Output and Continuous Dormand & Prince Pairs" pp. 188-191.
 *
*/

#define DIMENSIONS 2

namespace Dopri54
{
    // We will be using the Dormand-Prince 5(4) order method.
    // The Butcher Tableau is :
    constexpr double c2 = 1.0 / 5.0;
    constexpr double c3 = 3.0 / 10.0;
    constexpr double c4 = 4.0 / 5.0;
    constexpr double c5 = 8.0 / 9.0;
    constexpr double a21 = 1.0 / 5.0;
    constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    constexpr double a71 = 35.0 / 384.0, a72 = 500.0 / 1113.0, a73 = 125.0 / 192.0, a74 = -2187.0 / 6784.0, a75 = 11.0 / 84.0;

    constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    constexpr double e1 = 5179.0 / 57600.0, e3 = 7571.0 / 16695.0, e4 = 393.0 / 640.0, e5 = -92097.0 / 339200.0, e6 = 187.0 / 2100.0, e7 = 1.0 / 40.0;
}

namespace StepSizeController
{
    enum Controllers
    {
        STANDARD,
        H211PI,
        H312PID,
        H211B,
        H312B,
        PI42,

        /*
        ***Recommended Controllers with Stepsize Low-Pass Filters and their Problem Classes***
        *  Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
        *  https://dl.acm.org/doi/10.1145/641876.641877 -> Table III. page 24
        */

        /*--------------------------------------------------------------------------
         * kbeta1 | kbeta2 | kbeta3 | alpha2 | alpha3 | Class    | Problem Type
         *-------------------------------------------------------------------------
         * 1/b    | 1/b    | 0      | 1/b    | 0      | H211b    | medium to nonsmooth
         * 1/6    | 1/6    | 0      | 0      | 0      | H211 PI  | medium to nonsmooth
         * 1/b    | 2/b    | 1/b    | 3/b    | 1/b    | H312b    | nonsmooth
         * 1/18   | 1/9    | 1/18   | 0      | 0      | H312 PID | nonsmooth
         *-------------------------------------------------------------------------
         */
    };
}

class Solver
{

private:
    const double m_p = 4.0;         // the order corresponding to the RK method
    const double m_k = m_p + 1.0;   // EPS => k = p + 1 and EPUS => k = p
    const double m_kappa = 1.5;     // kappa ∈ [0.7, 2] as suggested in the literature
    const double m_acceptSF = 0.81; // accept safety factor

    std::function<void(double, std::array<double, DIMENSIONS> &, std::array<double, DIMENSIONS> &)> m_F{}; // f(t, y)

    std::array<double, DIMENSIONS> m_yn{}, m_X{}, m_K1{}, m_K2{}, m_K3{}, m_K4{}, m_K5{}, m_K6{}, m_K7{}, m_yn1{}, m_yn2{}, m_truncationErrors{}, m_sci{};

    double m_h{}, m_t{}, m_tFinal{}, m_absTol{}, m_relTol{};

    double m_beta1, m_beta2, m_beta3, m_alpha2, m_alpha3;

    bool m_denseOut;

    double m_cerr1 = 1.0, m_cerr2 = 1.0, m_rh1 = 1.0, m_rh2 = 1.0;

    int m_acceptedSteps{};
    int m_rejectedSteps{};

private:
    void SetControllerParameters(double b1, double b2, double b3, double a2, double a3)
    {
        m_beta1 = b1 / m_k;
        m_beta2 = b2 / m_k;
        m_beta3 = b3 / m_k;
        m_alpha2 = a2;
        m_alpha3 = a3;
    }

    double hairerNorm()
    {
        /**
         * ----- Calculate error-norm ||err|| -----
         * using the L2-Norm or the Euclidean Norm.
         * */
        double sum = 0.0;
        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            sum = sum + std::pow(m_truncationErrors[i] / m_sci[i], 2.0);
        }
        return std::sqrt(sum / static_cast<double>(DIMENSIONS));
    }

    double Process()
    {

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i];
        }

        m_F(m_t, m_X, m_K1); //--------------------- 1ST-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a21 * m_K1[i]);
        }

        m_F(m_t + Dopri54::c2 * m_h, m_X, m_K2); //--------------------- 2ND-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a31 * m_K1[i] + Dopri54::a32 * m_K2[i]);
        }

        m_F(m_t + Dopri54::c3 * m_h, m_X, m_K3); //--------------------- 3RD-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a41 * m_K1[i] + Dopri54::a42 * m_K2[i] +
                                      Dopri54::a43 * m_K3[i]);
        }

        m_F(m_t + Dopri54::c4 * m_h, m_X, m_K4); //--------------------- 4TH-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a51 * m_K1[i] + Dopri54::a52 * m_K2[i] +
                                      Dopri54::a53 * m_K3[i] + Dopri54::a54 * m_K4[i]);
        }

        m_F(m_t + Dopri54::c5 * m_h, m_X, m_K5); //--------------------- 5TH-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a61 * m_K1[i] + Dopri54::a62 * m_K2[i] +
                                      Dopri54::a63 * m_K3[i] + Dopri54::a64 * m_K4[i] + Dopri54::a65 * m_K5[i]);
        }

        m_F(m_t + m_h, m_X, m_K6); //--------------------- 6TH-stage -----------------------

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_X[i] = m_yn[i] + m_h * (Dopri54::a71 * m_K1[i] + Dopri54::a72 * m_K3[i] +
                                      Dopri54::a73 * m_K4[i] + Dopri54::a74 * m_K5[i] + Dopri54::a75 * m_K6[i]);
        }

        m_F(m_t + m_h, m_X, m_K7); //--------------------- 7TH-stage -----------------------

        // Calculate the 5th-order and 4th-order accurate solution.
        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_yn1[i] = m_yn[i] + m_h * (Dopri54::b1 * m_K1[i] + Dopri54::b3 * m_K3[i] +
                                        Dopri54::b4 * m_K4[i] + Dopri54::b5 * m_K5[i] + Dopri54::b6 * m_K6[i]); // 5th-Order accurate solution. Used to advance the solution.
            m_yn2[i] = m_yn[i] + m_h * (Dopri54::e1 * m_K1[i] + Dopri54::e3 * m_K3[i] +
                                        Dopri54::e4 * m_K4[i] + Dopri54::e5 * m_K5[i] + Dopri54::e6 * m_K6[i] + Dopri54::e7 * m_K7[i]); // 4th-Order accurate solution. Used for comparison to estimate error.
        }

        // Calculate local errors
        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            m_truncationErrors[i] = m_h * ((Dopri54::b1 - Dopri54::e1) * m_K1[i] + (Dopri54::b3 - Dopri54::e3) * m_K3[i] +
                                           (Dopri54::b4 - Dopri54::e4) * m_K4[i] + (Dopri54::b5 - Dopri54::e5) * m_K5[i] +
                                           (Dopri54::b6 - Dopri54::e6) * m_K6[i] - Dopri54::e7 * m_K7[i]);
        }

        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            // absTol and relTol are the desired tolerances prescribed by the user.
            m_sci[i] = m_absTol + std::max(std::fabs(m_yn[i]), std::fabs(m_yn1[i])) * m_relTol;
        }

        // local error is controlled by error-per-unit-steps (EPUS) or error-per-steps (EPS)
        // double r = L2Norm() / h;  // Error-per-unit-steps (EPUS)
        double err_estimate_scaled = hairerNorm(); // Error-per-steps (EPS)

        return err_estimate_scaled;
    }

    double Filter(double cerrPres, double cerrOld1, double cerrOld2, double rho1, double rho2)
    {
        /**
         * The General controller formula for Order-Dynamics pD <= 3 with Control error filtering.
         * Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
         * https://dl.acm.org/doi/10.1145/641876.641877 -> page 22
         */
        double result = std::pow(cerrPres, m_beta1) *
                        std::pow(cerrOld1, m_beta2) *
                        std::pow(cerrOld2, m_beta3) *
                        std::pow(rho1, -m_alpha2) *
                        std::pow(rho2, -m_alpha3);
        return result;
    }

    void ControlStepSize(double ratio)
    {
        m_h *= ratio;
    }

    std::array<double, DIMENSIONS> Interpolate(double theta, double hPresent)
    {

        const double C1 = 5.0 * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
        const double C3 = 100.0 * (882725551.0 - 15701508.0 * theta) / 32700410799.0;
        const double C4 = 25.0 * (443332067.0 - 31403016.0 * theta) / 1880347072.0;
        const double C5 = 32805.0 * (23143187.0 - 3489224.0 * theta) / 199316789632.0;
        const double C6 = 55.0 * (29972135.0 - 7076736.0 * theta) / 822651844.0;
        const double C7 = 10.0 * (7414447.0 - 829305.0 * theta) / 29380423.0;

        double theta_sqr = std::pow(theta, 2.0);
        double term1 = theta_sqr * (3.0 - 2.0 * theta);
        double term2 = theta_sqr * std::pow(theta - 1.0, 2.0);
        double term3 = theta * std::pow(theta - 1.0, 2.0);
        double term4 = (theta - 1.0) * std::pow(theta, 2.0);

        double b1Theta = term1 * Dopri54::b1 + term3 - term2 * C1;
        double b3Theta = term1 * Dopri54::b3 + term2 * C3;
        double b4Theta = term1 * Dopri54::b4 - term2 * C4;
        double b5Theta = term1 * Dopri54::b5 + term2 * C5;
        double b6Theta = term1 * Dopri54::b6 - term2 * C6;
        double b7Theta = term4 + term2 * C7;

        std::array<double, DIMENSIONS> solution{};
        for (size_t i = 0; i < DIMENSIONS; i++)
        {
            /* code */
            solution[i] = m_yn[i] + hPresent * (b1Theta * m_K1[i] + b3Theta * m_K3[i] + b4Theta * m_K4[i] +
                                                b5Theta * m_K5[i] + b6Theta * m_K6[i] + b7Theta * m_K7[i]);
        }
        return solution;
    }

public:
    std::vector<std::array<double, DIMENSIONS>> m_yOut{}; // accumulate solution steps
    std::vector<double> m_tOut{};                         // accumulate time steps

    // Solver(controllerType, f(t, y), y0, h, t0, tFinal, absTolerance=1E-6, relTolerance=1E-4, denseOut=false)
    Solver(StepSizeController::Controllers controller, std::function<void(double, std::array<double, DIMENSIONS> &, std::array<double, DIMENSIONS> &)> fName,
           std::array<double, DIMENSIONS> &y0, double stepSize, double t0, double tFinal, double absTol = 1e-6, double relTol = 1e-4, bool denseOut = false)
        : m_F{fName}, m_yn{y0}, m_h{stepSize}, m_t{t0}, m_tFinal{tFinal}, m_absTol{absTol}, m_relTol{relTol}, m_denseOut{denseOut}
    {

        switch (controller)
        {
        // set the parameters
        case StepSizeController::STANDARD:
            SetControllerParameters(1.0, 0.0, 0.0, 0.0, 0.0);
            break;
        case StepSizeController::H211PI:
            SetControllerParameters(1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0, 0.0);
            break;
        case StepSizeController::H312PID:
            SetControllerParameters(1.0 / 18.0, 1.0 / 9.0, 1.0 / 18.0, 0.0, 0.0);
            break;
        case StepSizeController::H211B: // b = 4.0
            SetControllerParameters(1.0 / 4.0, 1.0 / 4.0, 0.0, 1.0 / 4.0, 0.0);
            break;
        case StepSizeController::H312B: // b = 8.0
            SetControllerParameters(1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0);
            break;
        case StepSizeController::PI42:
            SetControllerParameters(3.0 / 5.0, -1.0 / 5.0, 0.0, 0.0, 0.0);
            break;
        }
    }

    void Solve()
    {
        // Add the initial value
        m_tOut.push_back(m_t);
        m_yOut.push_back(m_yn);

        // Closed-loop system
        while (m_t < m_tFinal)
        {
            m_h = std::min(m_h, m_tFinal - m_t);

            double err_estimate_scaled = Process();

            double cerr = 1.0 / err_estimate_scaled; // inverse of the error estimates scaled

            double rho = Filter(cerr, m_cerr1, m_cerr2, m_rh1, m_rh2);

            // Save previous values for the next step.
            m_cerr2 = m_cerr1;
            m_cerr1 = cerr;

            m_rh2 = m_rh1;
            m_rh1 = rho;

            // Apply a limiter
            double ratio = 1.0 + m_kappa * std::atan((rho - 1.0) / m_kappa);

            if (ratio < m_acceptSF)
            { // Reject steps and recalculate with the new stepsize
                ControlStepSize(ratio);
                m_rejectedSteps++;
                continue;
            }
            else
            { // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.

                if (m_denseOut) // do 1-more extra steps at the middle between yn and yn1.
                {
                    /*
                     * theta ∈ [0, 1], theta = 0 => yn, theta = 1 => yn1
                     * Interpolate at open-interval theta ∈ (0, 1)
                     * un+1(t + theta*h) = yn + h * sum(bi(theta)*Ki), i = 1...s, theta ∈ (0, 1)
                     */
                    double theta = 0.5;
                    std::array<double, DIMENSIONS> extraSteps = Interpolate(theta, m_h);
                    m_tOut.push_back(m_t + theta * m_h);
                    m_yOut.push_back(extraSteps);
                }

                m_t += m_h;
                ControlStepSize(ratio);
                m_acceptedSteps++;

                m_tOut.push_back(m_t);
                m_yOut.push_back(m_yn1);

                for (size_t i = 0; i < DIMENSIONS; i++)
                {
                    m_yn[i] = m_yn1[i];
                }
            }
        }
    }

    void DisplaySteps()
    {
        std::cout << "\n\tSteps: accepted = " << m_acceptedSteps << " rejected = " << m_rejectedSteps << std::endl;
    }

    void DisplayResults()
    {

        std::cout << std::endl;

        std::cout << std::setprecision(15) << std::fixed;

        for (size_t i = 0; i < m_yOut.size(); i++)
        {
            /* code */
            std::cout << "step " << i << " at t = " << m_tOut[i] << '\n';
            std::cout << "\t -> ";
            for (int j = 0; j < DIMENSIONS; j++)
            {
                std::cout << m_yOut[i][j] << ' ';
            }
            std::cout << '\n';
        }
    }
};

// Example problem: Van Der Pol Equation
// x' = y
// y' = mu * (1 - x*x) * y - x

// F(t, y)
void F(double t, std::array<double, DIMENSIONS> &X, std::array<double, DIMENSIONS> &Ks)
{
    double xDot = X[1];
    double yDot = 0.75 * (1.0 - X[0]*X[0]) * X[1] - X[0];
    Ks[0] = xDot;
    Ks[1] = yDot;
}

int main(void)
{

    std::array<double, DIMENSIONS> y0 = {2.0, 2.0};
    double t0 = 0.0;
    double tFinal = 10.0;
    double stepSize = 0.2;

    StepSizeController::Controllers controller = StepSizeController::H211PI; // Choose a controller
    Solver solver(controller, F, y0, stepSize, t0, tFinal, 1e-6, 1e-6);
    solver.Solve();
    solver.DisplayResults();
    solver.DisplaySteps();

    return EXIT_SUCCESS;
}