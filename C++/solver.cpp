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
 * Naval, Biliran, Philippines
 * 
 * References:
 * [1] K. Gustafsson. Control theoretic techniques for stepsize selection in explicit Runge–Kutta methods. ACM TOMS 17:533–554, 1991.
 * [2] G. Söderlind, “Digital Filters in Adaptive Time-Stepping”, ACM Trans. Math. Software, 29, pp. 1-26 (2003).
 * [3] Söderlind, Wang (2006) Adaptive time-stepping and computational stability DOI: 10.1016/j.cam.2005.03.008 
 * [4] Hairer, Ernst et. al. 𝘚𝘰𝘭𝘷𝘪𝘯𝘨 𝘖𝘳𝘥𝘪𝘯𝘢𝘳𝘺 𝘋𝘪𝘧𝘧𝘦𝘳𝘦𝘯𝘵𝘪𝘢𝘭 𝘌𝘲𝘶𝘢𝘵𝘪𝘰𝘯𝘴 𝘐: 𝘕𝘰𝘯𝘴𝘵𝘪𝘧𝘧 𝘗𝘳𝘰𝘣𝘭𝘦𝘮𝘴, Springer Science & Business Media, 1987, 
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


class Solver {

private :
    const double p         = 4.0;        // the order corresponding to the RK method
    const double k         = p + 1.0;    // EPS => k = p + 1 and EPUS => k = p
    const double kappa     = 1.35;       // kappa ∈ [0.7, 2] as suggested in the literature
    const double accept_SF = 0.90;       // accept safety factor

    double absTOL = 1e-6;    // default absolute tolerance
    double relTol = 1e-4;    // default relative tolerance

    double h, t, t_end;

    double beta1, beta2, beta3, alpha2, alpha3;

    bool dense_out = false;

    std::function<std::array<double, COMPONENTS>(double, std::array<double, COMPONENTS>&)> F{};   // f(t, y)

    std::array<double, COMPONENTS> yn{}, X{}, K1{}, K2{}, K3{}, K4{}, K5{}, K6{}, K7{}, yn1{}, yn2{}, local_errs{}, sci{};

    double cerr1 = 1.0, cerr2 = 1.0, rh1 = 1.0, rh2 = 1.0;

    int accepted_steps = 0;
    int rejected_steps = 0;

    std::vector<std::array<double, COMPONENTS>> y_out{};    // accumulate solution steps
    std::vector<double> t_out{};                            // accumulate time steps

private :
    void SetControllerParameters(double b1, double b2, double b3, double a2, double a3) {
        beta1  = b1/k;
        beta2  = b2/k;
        beta3  = b3/k;
        alpha2 = a2;
        alpha3 = a3;
    }

    double L2Norm(std::array<double, COMPONENTS> &err, std::array<double, COMPONENTS> &sci)
    {
        /**
         * ----- Calculate error-norm ||err|| -----
         * We will be using the L2-Norm or the Euclidean Norm.
         * */
        double sum = 0.0;
        for (size_t i = 0; i < COMPONENTS; i++)
        {
            sum = sum + pow(err[i] / sci[i], 2.0);
        }
        return sqrt(sum/COMPONENTS);
    }

    double ControlStepSize(double cerr_pres, double cerr_old1, double cerr_old2, double rho1, double rho2) {
        /**
         * The General controller formula for Order-Dynamics pD <= 3 with Control error filtering.
         * Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003) 
         * https://dl.acm.org/doi/10.1145/641876.641877 -> page 22
         */
        double result = pow(cerr_pres, beta1) * 
                        pow(cerr_old1, beta2) * 
                        pow(cerr_old2, beta3) * 
                        pow(rho1, -alpha2) * 
                        pow(rho2, -alpha3);
        return result;
    }

    double Interpolate(double theta, double h_present, double y0, double k1, double k3, double k4, double k5, double k6, double k7) {

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

        double b1_theta  = term1 * Dopri5::b1 + term3 - term2 * C1;
        double b3_theta  = term1 * Dopri5::b3 + term2 * C3;
        double b4_theta  = term1 * Dopri5::b4 - term2 * C4;
        double b5_theta  = term1 * Dopri5::b5 + term2 * C5;
        double b6_theta  = term1 * Dopri5::b6 - term2 * C6;
        double b7_theta  = term4 + term2 * C7;

        double solution = y0 + h_present * (b1_theta*k1 + b3_theta*k3 + b4_theta*k4 + b5_theta*k5 + b6_theta*k6 + b7_theta*k7);
        return solution;
    }

public :

    /* Setter Functions*/

    // Set 'true' to add extra steps to produce smoother steps.
    void SetDenseOut(bool is_dense_out) {
        dense_out = is_dense_out;
    }

    // Set the [A]bsolute and [R]elative tolerance
    void SetTolerance(double AbsTol, double RelTol) {
        absTOL = AbsTol;
        relTol = RelTol;
    }

    /*  Member Functions */
    void DisplayStatus() {
        std::cout << "Accepted steps: " << accepted_steps << " Rejected steps: " << rejected_steps << std::endl;
    }

    void DisplaySteps() {
        std::cout << std::setprecision(15) << std::fixed;

        for (size_t i = 0; i < y_out.size(); i++)
        {
            /* code */
            std::cout << "Step: " << i << " Time: " << t_out[i] << '\n';
            for (int j = 0; j < COMPONENTS; j++)
            {
                std::cout << y_out[i][j] << " ";
            }
            std::cout << "\n\n";
        }
    }

    Solver(int controller, std::array<double, COMPONENTS> &y0, double stepsize, double t0, double T_end, 
    std::function<std::array<double, COMPONENTS>(double, std::array<double, COMPONENTS>&)> fname)
    {
        t = t0;
        h = stepsize;
        t_end = T_end;
        F = fname;

        for (size_t i = 0; i < COMPONENTS; i++)
        {
            /* code */
            yn[i] = y0[i];
        }

        switch (controller) 
        {
            // set the parameters
            case 0 :   // STANDARD
                SetControllerParameters(1.0, 0.0, 0.0, 0.0, 0.0);
                break;
            case 1 :   // H211PI
                SetControllerParameters(1.0/6.0, 1.0/6.0, 0.0, 0.0, 0.0);
                break;
            case 2 :   // H312PID
                SetControllerParameters(1.0/18.0, 1.0/9.0, 1.0/18.0, 0.0, 0.0);
                break;
            case 3 :  // H211b, b = 4.0
                SetControllerParameters(1.0/4.0, 1.0/4.0, 0.0, 1.0/4.0, 0.0);
                break;
            case 4 :  //H312b, b = 8.0
                SetControllerParameters(1.0/8.0, 1.0/8.0, 1.0/8.0, 3.0/8.0, 1.0/8.0);
                break;
            case 5 :  // PI42
                SetControllerParameters(3.0/5.0, -1.0/5.0, 0.0, 0.0, 0.0);
                break;
        }
    }

    void Solve() {


        // Add the initial value
        t_out.push_back(t);
        y_out.push_back(yn);

        // Begin the stepper
        while (t < t_end) {
            h = std::min(h, t_end - t);

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i] = yn[i];
            }

            K1 = F(t, X); //--------------------- 1ST-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a21*K1[i]); 
            }

            K2 = F(t + Dopri5::c2*h, X); //--------------------- 2ND-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a31*K1[i] + Dopri5::a32*K2[i]);
            }

            K3 = F(t + Dopri5::c3*h, X); //--------------------- 3RD-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a41*K1[i] + Dopri5::a42*K2[i] + Dopri5::a43*K3[i]);
            }

            K4 = F(t + Dopri5::c4*h, X); //--------------------- 4TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a51*K1[i] + Dopri5::a52*K2[i] + Dopri5::a53*K3[i] + Dopri5::a54*K4[i]);
            }

            K5 = F(t + Dopri5::c5*h, X); //--------------------- 5TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a61*K1[i] + Dopri5::a62*K2[i] + Dopri5::a63*K3[i] + Dopri5::a64*K4[i] + Dopri5::a65*K5[i]);
            }

            K6 = F(t + h, X); //--------------------- 6TH-stage -----------------------

            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                X[i]  = yn[i] + h * (Dopri5::a71*K1[i] + Dopri5::a72*K3[i] + Dopri5::a73*K4[i] + Dopri5::a74*K5[i] + Dopri5::a75*K6[i]);
            }

            K7 = F(t + h, X); //--------------------- 7TH-stage -----------------------

            // Calculate the 5th-order accurate solution and the alternative 4th-order solution for the error estimation.
            for (size_t i = 0; i < COMPONENTS; i++)
            {
                yn1[i] = yn[i] + h * (Dopri5::b1*K1[i] + Dopri5::b3*K3[i] + Dopri5::b4*K4[i] + Dopri5::b5*K5[i] + Dopri5::b6*K6[i]);               // 5th-Order accurate solution.
                yn2[i] = yn[i] + h * (Dopri5::e1*K1[i] + Dopri5::e3*K3[i] + Dopri5::e4*K4[i] + Dopri5::e5*K5[i] + Dopri5::e6*K6[i] + Dopri5::e7*K7[i]);   // 4th-Order alternative solution for local-error estimation.
            }

            // Calculate local errors
            for (size_t i = 0; i < COMPONENTS; i++) 
            {
                local_errs[i] = h * (yn2[i] - yn1[i]);
            }

            for (size_t i = 0; i < COMPONENTS; i++)
            {
                /* 
                 * absTol and relTol are the desired tolerances prescribed by the user. 
                 * Note: absTol and relTol can be a vector meaning you can set a different values of absTol and relTol for every component in yn.
                 * Ex. If we have a 2D(x, y)-component vector then x = (absTol, relTol) = (1e-06, 1e-03) and y = (absTol, relTol) = (1e-08, 1e-06).
                 */
                sci[i] = absTOL + std::max(fabs(yn[i]), fabs(yn1[i])) * relTol;  // We will just use the same absTol and relTol for all component in yn.
            }

            // local error is controlled by error-per-unit-steps (EPUS) or error-per-steps (EPS)
            // let r = L2Norm(local_errs, sci) / h;  // Error-per-unit-steps (EPUS)
            double r = L2Norm(local_errs, sci);   // Error-per-steps (EPS)

            double cerr = 1.0 / r;

            double rho = ControlStepSize(cerr, cerr1, cerr2, rh1, rh2);

            // Save previous values for the next step.
            cerr2 = cerr1;
            cerr1 = cerr;

            rh2 = rh1;
            rh1 = rho;

            // Apply a limiter
            double ratio = 1.0 + kappa * atan((rho - 1.0) / kappa);

            if (ratio < accept_SF) { // Reject steps and recalculate with the new stepsize
                h = ratio * h;
                rejected_steps++;
                continue;
            } 
            else {  // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.

                if (dense_out) {
                    
                    std::array<double, COMPONENTS> extra_steps{};
                    double theta = 0.5; /// theta = [0, 1], theta = 0 => yn, theta = 1 => yn1
                    for (int i = 0; i < COMPONENTS; i++)
                    {
                        // interpolate at theta = 0.5
                        double step    = Interpolate(0.5, h, yn[i], K1[i], K3[i], K4[i], K5[i], K6[i], K7[i]);     // Interpolate on the ith-component of yn
                        extra_steps[i] = step;
                    }
                    t_out.push_back(t + theta*h);
                    y_out.push_back(extra_steps);
                }

                t = t + h;
                h = ratio * h;
                accepted_steps++;

                t_out.push_back(t);
                y_out.push_back(yn1);

                for (size_t i = 0; i < COMPONENTS; i++)
                {
                    yn[i] = yn1[i];
                }
            }
        }
    }
};

enum class Controller
{
    // Warning: Don't change the order.
    STANDARD, // assign 0
    H211PI,   // assign 1
    H312PID,  // assign 2
    H211B,    // assign 3
    H312B,    // assign 4
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
    double stepsize = 0.2;
    double T_end = 2.0;

    Controller controller = Controller::H312B;      // Choose a controller

    Solver solver(static_cast<int> (controller), y0, stepsize, t0, T_end, F);
    solver.SetTolerance(1e-06, 1e-04);
    solver.SetDenseOut(false);
    solver.Solve();
    solver.DisplaySteps();
    solver.DisplayStatus();
    return EXIT_SUCCESS;
}