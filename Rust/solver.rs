/*
 * Copyright (c) 2024 Mckenzie Regalado - All Rights Reserved.
 * Naval, Biliran, Philippines
 *
 * References:
 * [1] K. Gustafsson. Control theoretic techniques for stepsize selection in explicit Runge–Kutta methods. ACM TOMS 17:533–554, 1991.
 * [2] G. Söderlind, “Digital Filters in Adaptive Time-Stepping”, ACM Trans. Math. Software, 29, pp. 1-26 (2003).
 * [3] Söderlind, Wang (2006) Adaptive time-stepping and computational stability DOI: 10.1016/j.cam.2005.03.008
 * [4] Hairer, Ernst et. al. 𝘚𝘰𝘭𝘷𝘪𝘯𝘨 𝘖𝘳𝘥𝘪𝘯𝘢𝘳𝘺 𝘋𝘪𝘧𝘧𝘦𝘳𝘦𝘯𝘵𝘪𝘢𝘭 𝘌𝘲𝘶𝘢𝘵𝘪𝘰𝘯𝘴 𝘐: 𝘕𝘰𝘯𝘴𝘵𝘪𝘧𝘧 𝘗𝘳𝘰𝘣𝘭𝘦𝘮𝘴, Springer Science & Business Media, 1987,
 * "Automatic Step Size Control" page-167, "Dense Output and Continuous Dormand & Prince Pairs" pp. 188-191.
*/

const COMPONENTS: usize = 1;

// Dormand-Prince 5(4) method butcher tableau
const C2: f64 = 1.0 / 5.0;
const C3: f64 = 3.0 / 10.0;
const C4: f64 = 4.0 / 5.0;
const C5: f64 = 8.0 / 9.0;
const A21: f64 = 1.0 / 5.0;
const A31: f64 = 3.0 / 40.0;
const A32: f64 = 9.0 / 40.0;
const A41: f64 = 44.0 / 45.0;
const A42: f64 = -56.0 / 15.0;
const A43: f64 = 32.0 / 9.0;
const A51: f64 = 19372.0 / 6561.0;
const A52: f64 = -25360.0 / 2187.0;
const A53: f64 = 64448.0 / 6561.0;
const A54: f64 = -212.0 / 729.0;
const A61: f64 = 9017.0 / 3168.0;
const A62: f64 = -355.0 / 33.0;
const A63: f64 = 46732.0 / 5247.0;
const A64: f64 = 49.0 / 176.0;
const A65: f64 = -5103.0 / 18656.0;
const A71: f64 = 35.0 / 384.0;
const A72: f64 = 500.0 / 1113.0;
const A73: f64 = 125.0 / 192.0;
const A74: f64 = -2187.0 / 6784.0;
const A75: f64 = 11.0 / 84.0;

const B1: f64 = 35.0 / 384.0;
const B3: f64 = 500.0 / 1113.0;
const B4: f64 = 125.0 / 192.0;
const B5: f64 = -2187.0 / 6784.0;
const B6: f64 = 11.0 / 84.0;
const E1: f64 = 5179.0 / 57600.0;
const E3: f64 = 7571.0 / 16695.0;
const E4: f64 = 393.0 / 640.0;
const E5: f64 = -92097.0 / 339200.0;
const E6: f64 = 187.0 / 2100.0;
const E7: f64 = 1.0 / 40.0;

// Set the parameters for controllers.

/*
 Recommended Controllers with Stepsize Low-Pass Filters and their Problem Classes
 Reference: Digital Filters in Adaptive Time-Stepping (Sorderlind, 2003)
 https://dl.acm.org/doi/10.1145/641876.641877 -> Table III. page 24

 *--------------------------------------------------------------------------
 * kbeta1 | kbeta2 | kbeta3 | alpha2 | alpha3 | Class    | Problem Type
 *-------------------------------------------------------------------------
 * 1/b    | 1/b    | 0      | 1/b    | 0      | H211b    | medium to nonsmooth
 * 1/6    | 1/6    | 0      | 0      | 0      | H211 PI  | medium to nonsmooth
 * 1/b    | 2/b    | 1/b    | 3/b    | 1/b    | H312b    | nonsmooth
 * 1/18   | 1/9    | 1/18   | 0      | 0      | H312 PID | nonsmooth
 *-------------------------------------------------------------------------
*/

/* For Standard controller */
const BETA1: f64 = 1.0;
const BETA2: f64 = 0.0;
const BETA3: f64 = 0.0;
const ALPHA2: f64 = 0.0;
const ALPHA3: f64 = 0.0;

/* For H211PI controller */
// const BETA1: f64 = 1.0 / 6.0;
// const BETA2: f64 = 1.0 / 6.0;
// const BETA3: f64 = 0.0;
// const ALPHA2: f64 = 0.0;
// const ALPHA3: f64 = 0.0;

/* For H312PID controller */
// const BETA1: f64 = 1.0 / 18.0;
// const BETA2: f64 = 1.0 / 9.0;
// const BETA3: f64 = 1.0 / 18.0;
// const ALPHA2: f64 = 0.0;
// const ALPHA3: f64 = 0.0;

/* For H211B controller, B = 4.0 */
// const BETA1: f64 = 1.0 / 4.0;
// const BETA2: f64 = 1.0 / 4.0;
// const BETA3: f64 = 0.0;
// const ALPHA2: f64 = 1.0 / 4.0;
// const ALPHA3: f64 = 0.0;

/* For H312B controller, B = 8.0 */
// const BETA1: f64 = 1.0 / 8.0;
// const BETA2: f64 = 1.0 / 8.0;
// const BETA3: f64 = 1.0 / 8.0;
// const ALPHA2: f64 = 3.0 / 8.0;
// const ALPHA3: f64 = 1.0 / 8.0;

/* For PI42 controller */
// const BETA1: f64 = 3.0 / 5.0;
// const BETA2: f64 = -1.0 / 5.0;
// const BETA3: f64 = 0.0;
// const ALPHA2: f64 = 0.0;
// const ALPHA3: f64 = 0.0;

const K: f64 = 4.0; // EPS => k = p + 1 and EPUS => k = p, where p is the order-of-error of the method.
const KAPPA: f64 = 1.35; // kappa ∈ [0.7, 2] as suggested in the literature
const ACCEPT_SF: f64 = 0.90; // accept safety factor

// ODE function
// Example Problem
// y' = y - t^2 + 1
// y(0) = 0.5
// exact solution: t**2 + 2t + 1 - 0.5e**t
fn f(t: f64, y: &[f64; COMPONENTS]) -> [f64; COMPONENTS] {
    // y' = y - t^2 + 1
    let mut dydt = [0.0; COMPONENTS];
    dydt[0] = y[0] - t * t + 1.0;
    return dydt;
}

fn l2_norm(err: &[f64; COMPONENTS], sc_i: &[f64; COMPONENTS]) -> f64 {
    /*
     * ----- Calculate error-norm ||err|| -----
     * We will be using the L2-Norm or the Euclidean Norm.
     *
     */
    let mut sum: f64 = 0.0;
    for i in 0..COMPONENTS {
        sum = sum + (err[i] / sc_i[i]).powf(2.0);
    }
    let n: f64 = COMPONENTS as f64;
    return (sum / n).sqrt();
}

fn control_step_size(cerr_press: f64, cerr_old1: f64, cerr_old2: f64, rho1: f64, rho2: f64) -> f64 {
    let result: f64 = cerr_press.powf(BETA1 / K)
        * cerr_old1.powf(BETA2 / K)
        * cerr_old2.powf(BETA3 / K)
        * rho1.powf(-ALPHA2)
        * rho2.powf(-ALPHA3);
    return result;
}

fn interpolate(
    theta: f64,
    h_present: f64,
    y_present: f64,
    k1: f64,
    k3: f64,
    k4: f64,
    k5: f64,
    k6: f64,
    k7: f64,
) -> f64 {
    let c1: f64 = 5.0 * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
    let c3: f64 = 100.0 * (882725551.0 - 15701508.0 * theta) / 32700410799.0;
    let c4: f64 = 25.0 * (443332067.0 - 31403016.0 * theta) / 1880347072.0;
    let c5: f64 = 32805.0 * (23143187.0 - 3489224.0 * theta) / 199316789632.0;
    let c6: f64 = 55.0 * (29972135.0 - 7076736.0 * theta) / 822651844.0;
    let c7: f64 = 10.0 * (7414447.0 - 829305.0 * theta) / 29380423.0;

    let theta_sqr = theta.powf(2.0);
    let term1 = theta_sqr * (3.0 - 2.0 * theta);
    let term2 = theta_sqr * (theta - 1.0).powf(2.0);
    let term3 = theta * (theta - 1.0).powf(2.0);
    let term4 = (theta - 1.0) * theta.powf(2.0);

    let b1_theta = term1 * B1 + term3 - term2 * c1;
    let b3_theta = term1 * B3 + term2 * c3;
    let b4_theta = term1 * B4 - term2 * c4;
    let b5_theta = term1 * B5 + term2 * c5;
    let b6_theta = term1 * B6 - term2 * c6;
    let b7_theta = term4 + term2 * c7;

    let solution = y_present
        + h_present
            * (b1_theta * k1
                + b3_theta * k3
                + b4_theta * k4
                + b5_theta * k5
                + b6_theta * k6
                + b7_theta * k7);
    return solution;
}

fn solve(
    t0: f64,
    y0: [f64; COMPONENTS],
    stepsize: f64,
    tfinal: f64,
    tol: (f64, f64),
    dense_out: bool,
    y_out: &mut Vec<[f64; COMPONENTS]>,
    t_out: &mut Vec<f64>,
) {
    let mut t = t0;
    let mut y = y0;
    let mut h = stepsize;
    let (abs_tol, rel_tol) = tol;

    let mut k1: [f64; COMPONENTS];
    let mut k2: [f64; COMPONENTS];
    let mut k3: [f64; COMPONENTS];
    let mut k4: [f64; COMPONENTS];
    let mut k5: [f64; COMPONENTS];
    let mut k6: [f64; COMPONENTS];
    let mut k7: [f64; COMPONENTS];
    let mut x: [f64; COMPONENTS] = [0.0; COMPONENTS];
    let mut y1: [f64; COMPONENTS] = [0.0; COMPONENTS];
    let mut y2: [f64; COMPONENTS] = [0.0; COMPONENTS];
    let mut local_errs: [f64; COMPONENTS] = [0.0; COMPONENTS];
    let mut sc_i: [f64; COMPONENTS] = [0.0; COMPONENTS];

    let mut cerr1: f64 = 1.0;
    let mut cerr2: f64 = 1.0;
    let mut rh1: f64 = 1.0;
    let mut rh2: f64 = 1.0;

    let mut accepted_steps: i32 = 0;
    let mut rejected_steps: i32 = 0;

    // Main stepper
    while t < tfinal {
        h = f64::min(h, tfinal - t);

        for i in 0..COMPONENTS {
            x[i] = y[i];
        }

        // 1st-stage
        k1 = f(t, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A21 * k1[i]);
        }

        // 2nd-stage
        k2 = f(t + C2 * h, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
        }

        // 3rd-stage
        k3 = f(t + C3 * h, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
        }

        // 4th-stage
        k4 = f(t + C4 * h, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
        }

        // 5th-stage
        k5 = f(t + C5 * h, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
        }

        // 6th-stage
        k6 = f(t + h, &x);
        for i in 0..COMPONENTS {
            x[i] = y[i] + h * (A71 * k1[i] + A72 * k3[i] + A73 * k4[i] + A74 * k5[i] + A75 * k6[i]);
        }

        k7 = f(t + h, &x);
        for i in 0..COMPONENTS {
            // 5th-Order accurate solution.
            y1[i] = y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]);
            // 4th-Order alternative solution for local-error estimation.
            y2[i] = y[i]
                + h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
        }

        // calculate local errors.
        for i in 0..COMPONENTS {
            local_errs[i] = h * (y2[i] - y1[i]);
        }

        for i in 0..COMPONENTS {
            /*
             * absTol and relTol are the desired tolerances prescribed by the user.
             * Note: absTol and relTol can be a vector meaning you can set a different values of absTol and relTol for every component in y.
             * Ex. If we have a 2D(x, y)-component vector then we can set x = (absTol, relTol) = (1e-06, 1e-03) and
             * y = (absTol, relTol) = (1e-08, 1e-06).
             */
            sc_i[i] = abs_tol + f64::max(y[i].abs(), y1[i].abs()) * rel_tol; // We will just use the same absTol and relTol for all component in y.
        }
        let err_norm: f64 = l2_norm(&local_errs, &sc_i);
        let r: f64 = err_norm / h; // local error is controlled by error per unit step (EPUS).
        let cerr: f64 = 1.0 / r;
        let rho: f64 = control_step_size(cerr, cerr1, cerr2, rh1, rh2);

        // Save previous values for the next step.
        cerr2 = cerr1;
        cerr1 = cerr;
        rh2 = rh1;
        rh1 = rho;

        let ratio: f64 = 1.0 + KAPPA * f64::atan((rho - 1.0) / KAPPA);

        // Reject steps and recalculate with the new stepsize
        if ratio < ACCEPT_SF {
            h = ratio * h;
            rejected_steps += 1;
            continue;
        }
        // Accept steps and the solution is advanced with yn1 and tried with the new stepsize.
        else {
            if dense_out {
                let mut extra_steps: [f64; COMPONENTS] = [0.0; COMPONENTS];
                let theta: f64 = 0.5; // theta = [0, 1], theta = 0 => y, theta = 1 => y1
                for i in 0..COMPONENTS {
                    // interpolate at theta = 0.5
                    extra_steps[i] =
                        interpolate(theta, h, y[i], k1[i], k3[i], k4[i], k5[i], k6[i], k7[i]);
                }
                y_out.push(extra_steps);
                t_out.push(t + theta * h);
            }
            t = t + h;
            h = ratio * h;
            accepted_steps += 1;
            for i in 0..COMPONENTS {
                y[i] = y1[i];
            }
            y_out.push(y);
            t_out.push(t);
        }
    }
    println!(
        "Accepted steps: {accepted_steps} Rejected steps: {rejected_steps} dense out? {dense_out}"
    );
}

fn display_steps(y_out: &Vec<[f64; COMPONENTS]>, t_out: &Vec<f64>) {
    for i in 0..y_out.len() {
        println!("Time: {}", t_out[i]);
        for component in 0..COMPONENTS {
            println!("Component {}: {}", component, y_out[i][component]);
        }
    }
}

fn main() {
    let y0: [f64; COMPONENTS] = [0.5];
    let t0: f64 = 0.0;

    let mut y_out: Vec<[f64; COMPONENTS]> = Vec::new(); // accumulate solution steps.
    let mut t_out: Vec<f64> = Vec::new(); // accumulate time steps.

    y_out.push(y0);
    t_out.push(t0);

    solve(
        t0,
        y0,
        0.2,
        2.0,
        (1e-06, 1e-04),
        false,
        &mut y_out,
        &mut t_out,
    );
    display_steps(&y_out, &t_out);
}
