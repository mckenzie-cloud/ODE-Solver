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

#define COMPONENTS 1 // It is simply the size of the array (ex. 2D-array is (x, y) component, 3D-array is (x, y, z) component

namespace Dopri5
{
    // We will be using the Dormand-Prince 5(4) order method.
    // The Butcher Tableau is :
    constexpr double c2  = 1.0/5.0;
    constexpr double c3  = 3.0/10.0;
    constexpr double c4  = 4.0/5.0;
    constexpr double c5  = 8.0/9.0;
    constexpr double a21 = 1.0/5.0;
    constexpr double a31 = 3.0/40.0,       a32 = 9.0/40.0;
    constexpr double a41 = 44.0/45.0,      a42 = -56.0/15.0,      a43 = 32.0/9.0;
    constexpr double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    constexpr double a61 = 9017.0/3168.0,  a62 = -355.0/33.0,     a63 = 46732.0/5247.0, a64 = 49.0/176.0,        a65 = -5103.0/18656.0;
    constexpr double a71 = 35.0/384.0,     a72 = 500.0/1113.0,    a73 = 125.0/192.0,    a74 = -2187.0/6784.0,    a75 = 11.0/84.0;

    constexpr double b1  = 35.0/384.0,      b3 = 500.0/1113.0,    b4  = 125.0/192.0,    b5  = -2187.0/6784.0,    b6  = 11.0/84.0;
    constexpr double e1  = 5179.0/57600.0,  e3 = 7571.0/16695.0,  e4  = 393.0/640.0,    e5  = -92097.0/339200.0, e6  = 187.0/2100.0, e7 = 1.0/40.0;
}

namespace Controller
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


class Solver {

private :
    const double m_p         = 4.0;        // the order corresponding to the RK method
    const double m_k         = m_p + 1.0;    // EPS => k = p + 1 and EPUS => k = p
    const double m_kappa     = 1.0;        // kappa ∈ [0.7, 2] as suggested in the literature
    const double m_acceptSF = 0.90;       // accept safety factor

    double m_h{}, m_t{}, m_tFinal{}, m_absTol{}, m_relTol{};

    double m_beta1, m_beta2, m_beta3, m_alpha2, m_alpha3;

    bool m_denseOut;

    std::function<std::array<double, COMPONENTS>(double, std::array<double, COMPONENTS>&)> m_F{};   // f(t, y)

    std::array<double, COMPONENTS> m_yn{}, m_X{}, m_K1{}, m_K2{}, m_K3{}, m_K4{}, m_K5{}, m_K6{}, m_K7{}, m_yn1{}, m_yn2{}, m_localErrs{}, m_sci{};

    double m_cerr1 = 1.0, m_cerr2 = 1.0, m_rh1 = 1.0, m_rh2 = 1.0;

    int m_acceptedSteps{};
    int m_rejectedSteps{};

    std::vector<std::array<double, COMPONENTS>> m_yOut{};    // accumulate solution steps
    std::vector<double> m_tOut{};                            // accumulate time steps

private :
    void SetControllerParameters(double b1, double b2, double b3, double a2, double a3) {
        m_beta1  = b1/m_k;
        m_beta2  = b2/m_k;
        m_beta3  = b3/m_k;
        m_alpha2 = a2;
        m_alpha3 = a3;
    }

    double L2Norm()
    {
        /**
         * ----- Calculate error-norm ||err|| -----
         * We will be using the L2-Norm or the Euclidean Norm.
         * */
        double sum = 0.0;
        for (size_t i = 0; i < COMPONENTS; i++)
        {
            sum = sum + pow(m_localErrs[i] / m_sci[i], 2.0);
        }
        return sqrt(sum / static_cast<double>(COMPONENTS));
    }

    double ControlStepSize(double cerrPres, double cerrOld1, double cerrOld2, double rho1, double rho2) {
        /**
         * The General controller formula for Order-Dynamics pD <= 3 with Control error filtering.
         * Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003) 
         * https://dl.acm.org/doi/10.1145/641876.641877 -> page 22
         */
        double result = pow(cerrPres, m_beta1) * 
                        pow(cerrOld1, m_beta2) * 
                        pow(cerrOld2, m_beta3) * 
                        pow(rho1, -m_alpha2) * 
                        pow(rho2, -m_alpha3);
        return result;
    }

    std::array<double, COMPONENTS> Interpolate(double theta, double hPresent) {

        const double C1  = 5.0     * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
        const double C3  = 100.0   * (882725551.0  - 15701508.0 * theta) / 32700410799.0;
        const double C4  = 25.0    * (443332067.0  - 31403016.0 * theta) / 1880347072.0;
        const double C5  = 32805.0 * (23143187.0   - 3489224.0  * theta) / 199316789632.0;
        const double C6  = 55.0    * (29972135.0   - 7076736.0  * theta) / 822651844.0;
        const double C7  = 10.0    * (7414447.0    - 829305.0   * theta) / 29380423.0;

        double theta_sqr = pow(theta, 2.0);
        double term1     = theta_sqr * (3.0 - 2.0 * theta);
        double term2     = theta_sqr * pow(theta - 1.0, 2.0);
        double term3     = theta     * pow(theta - 1.0, 2.0);
        double term4     = (theta - 1.0) * pow(theta, 2.0);

        double b1Theta  = term1 * Dopri5::b1 + term3 - term2 * C1;
        double b3Theta  = term1 * Dopri5::b3 + term2 * C3;
        double b4Theta  = term1 * Dopri5::b4 - term2 * C4;
        double b5Theta  = term1 * Dopri5::b5 + term2 * C5;
        double b6Theta  = term1 * Dopri5::b6 - term2 * C6;
        double b7Theta  = term4 + term2 * C7;

        std::array<double, COMPONENTS> solution{};
        for (size_t i = 0; i < COMPONENTS; i++)
        {
            /* code */
            solution[i] = m_yn[i] + hPresent * (b1Theta*m_K1[i] + b3Theta*m_K3[i] + b4Theta*m_K4[i] + 
                                                b5Theta*m_K5[i] + b6Theta*m_K6[i] + b7Theta*m_K7[i]);
        }
        return solution;
    }

public :

    /* Setter Functions*/
    /*  Member Functions */
    void DisplayStatus() {
        std::cout << "Accepted steps: " << m_acceptedSteps << " Rejected steps: " << m_rejectedSteps << std::endl;
    }

    void DisplaySteps() {
        std::cout << std::setprecision(15) << std::fixed;

        for (size_t i = 0; i < m_yOut.size(); i++)
        {
            /* code */
            std::cout << "Step: " << i << " Time: " << m_tOut[i] << '\n';
            for (int j = 0; j < COMPONENTS; j++)
            {
                std::cout << m_yOut[i][j] << " ";
            }
            std::cout << "\n\n";
        }
    }

    // Solver(controllerType, f(t, y), y0, h, t0, tFinal, absTolerance=1E-6, relTolerance=1E-4, denseOut=false)
    Solver(Controller::Controllers controller, std::function<std::array<double, COMPONENTS>(double, std::array<double, COMPONENTS>&)> fName, 
           std::array<double, COMPONENTS> &y0, double stepSize, double t0, double tFinal, double absTol=1e-6, double relTol=1e-4, bool denseOut=false) : 
           m_F { fName },
           m_t { t0 },  
           m_h { stepSize }, 
           m_tFinal { tFinal }, 
           m_absTol { absTol }, 
           m_relTol { relTol },
           m_denseOut { denseOut }
    {

        for (size_t i = 0; i < COMPONENTS; i++)
        {
            /* code */
            m_yn[i] = y0[i];
        }

        switch (controller) 
        {
            // set the parameters
            case Controller::STANDARD :   
                SetControllerParameters(1.0, 0.0, 0.0, 0.0, 0.0);
                break;
            case Controller::H211PI :   
                SetControllerParameters(1.0/6.0, 1.0/6.0, 0.0, 0.0, 0.0);
                break;
            case Controller::H312PID :   
                SetControllerParameters(1.0/18.0, 1.0/9.0, 1.0/18.0, 0.0, 0.0);
                break;
            case Controller::H211B :  // b = 4.0
                SetControllerParameters(1.0/4.0, 1.0/4.0, 0.0, 1.0/4.0, 0.0);
                break;
            case Controller::H312B :  // b = 8.0
                SetControllerParameters(1.0/8.0, 1.0/8.0, 1.0/8.0, 3.0/8.0, 1.0/8.0);
                break;
            case Controller::PI42 :
                SetControllerParameters(3.0/5.0, -1.0/5.0, 0.0, 0.0, 0.0);
                break;
        }
    }

    void Solve() 
    {
        // Add the initial value
        m_tOut.push_back(m_t);
        m_yOut.push_back(m_yn);

        // Begin the stepper
        while (m_t < m_tFinal) 
        {
            m_h = std::min(m_h, m_tFinal - m_t);

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i];
            }

            m_K1 = m_F(m_t, m_X); //--------------------- 1ST-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a21*m_K1[i]); 
            }

            m_K2 = m_F(m_t + Dopri5::c2*m_h, m_X); //--------------------- 2ND-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a31*m_K1[i] + Dopri5::a32*m_K2[i]);
            }

            m_K3 = m_F(m_t + Dopri5::c3*m_h, m_X); //--------------------- 3RD-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a41*m_K1[i] + Dopri5::a42*m_K2[i] + 
                                          Dopri5::a43*m_K3[i]);
            }

            m_K4 = m_F(m_t + Dopri5::c4*m_h, m_X); //--------------------- 4TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a51*m_K1[i] + Dopri5::a52*m_K2[i] + 
                                          Dopri5::a53*m_K3[i] + Dopri5::a54*m_K4[i]);
            }

            m_K5 = m_F(m_t + Dopri5::c5*m_h, m_X); //--------------------- 5TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a61*m_K1[i] + Dopri5::a62*m_K2[i] + 
                                          Dopri5::a63*m_K3[i] + Dopri5::a64*m_K4[i] + Dopri5::a65*m_K5[i]);
            }

            m_K6 = m_F(m_t + m_h, m_X); //--------------------- 6TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_X[i] = m_yn[i] + m_h * (Dopri5::a71*m_K1[i] + Dopri5::a72*m_K3[i] + 
                                          Dopri5::a73*m_K4[i] + Dopri5::a74*m_K5[i] + Dopri5::a75*m_K6[i]);
            }

            m_K7 = m_F(m_t + m_h, m_X); //--------------------- 7TH-stage -----------------------

            // Calculate the 5th-order and 4th-order accurate solution.
            for (size_t i = 0; i < COMPONENTS; i++)
            {
                m_yn1[i] = m_yn[i] + m_h * (Dopri5::b1*m_K1[i] + Dopri5::b3*m_K3[i] + 
                                            Dopri5::b4*m_K4[i] + Dopri5::b5*m_K5[i] + Dopri5::b6*m_K6[i]);     // 5th-Order accurate solution. Used to advance the solution.
                m_yn2[i] = m_yn[i] + m_h * (Dopri5::e1*m_K1[i] + Dopri5::e3*m_K3[i] + 
                                            Dopri5::e4*m_K4[i] + Dopri5::e5*m_K5[i] + Dopri5::e6*m_K6[i] + Dopri5::e7*m_K7[i]);   // 4th-Order accurate solution. Used for comparison to estimate error.
            }

            // Calculate local errors
            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                m_localErrs[i] = m_h * ((Dopri5::b1 - Dopri5::e1)*m_K1[i] + 
                                        (Dopri5::b3 - Dopri5::e3)*m_K3[i] +
                                        (Dopri5::b4 - Dopri5::e4)*m_K4[i] +
                                        (Dopri5::b5 - Dopri5::e5)*m_K5[i] +
                                        (Dopri5::b6 - Dopri5::e6)*m_K6[i] - Dopri5::e7*m_K7[i]);
            }

            for (size_t i = 0; i < COMPONENTS; i++)
            {
                // absTol and relTol are the desired tolerances prescribed by the user. 
                m_sci[i] = m_absTol + std::max(fabs(m_yn[i]), fabs(m_yn1[i])) * m_relTol;
            }

            // local error is controlled by error-per-unit-steps (EPUS) or error-per-steps (EPS)
            // let r = L2Norm(localErrs, sci) / h;  // Error-per-unit-steps (EPUS)
            double r = L2Norm();   // Error-per-steps (EPS)

            double cerr = 1.0 / r;

            double rho = ControlStepSize(cerr, m_cerr1, m_cerr2, m_rh1, m_rh2);

            // Save previous values for the next step.
            m_cerr2 = m_cerr1;
            m_cerr1 = cerr;

            m_rh2 = m_rh1;
            m_rh1 = rho;

            // Apply a limiter
            double ratio = 1.0 + m_kappa * atan((rho - 1.0) / m_kappa);

            if (ratio < m_acceptSF) { // Reject steps and recalculate with the new stepsize
                m_h = ratio * m_h;
                m_rejectedSteps++;
                continue;
            } 
            else {  // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.

                if (m_denseOut) 
                {
                    /*
                     * theta ∈ [0, 1], theta = 0 => yn, theta = 1 => yn1
                     * Interpolate at open-interval theta ∈ (0, 1)
                     * un+1(t + theta*h) = yn + h * sum(bi(theta)*Ki), i = 1...s, theta ∈ (0, 1)
                     */
                    double theta = 0.5; 
                    std::array<double, COMPONENTS> extraSteps = Interpolate(theta, m_h);
                    m_tOut.push_back(m_t + theta*m_h);
                    m_yOut.push_back(extraSteps);
                }

                m_t = m_t + m_h;
                m_h = ratio * m_h;
                m_acceptedSteps++;

                m_tOut.push_back(m_t);
                m_yOut.push_back(m_yn1);

                for (size_t i = 0; i < COMPONENTS; i++)
                {
                    m_yn[i] = m_yn1[i];
                }
            }
        }
    }
};

// Example problem
// y' = y - t^2 + 1
// y(0) = 0.5
// exact solution: t**2 + 2t + 1 - 0.5e**t
std::array<double, COMPONENTS> F(double t, std::array<double, COMPONENTS> &X) {    // F(t, y)
    double dotx = X[0] - t*t + 1.0;
    return std::array<double, COMPONENTS> {dotx};
}

int main(void) {

    std::array<double, COMPONENTS> y0 = {0.5};
    double t0 = 0.0;
    double tFinal = 2.0;
    double stepSize = 0.2;

    Controller::Controllers controller = Controller::STANDARD;      // Choose a controller
    Solver solver(controller, F, y0, stepSize, t0, tFinal);
    solver.Solve();
    solver.DisplaySteps();
    solver.DisplayStatus();
    return EXIT_SUCCESS;
}