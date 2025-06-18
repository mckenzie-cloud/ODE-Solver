
namespace Dopri54
{
    // We will be using the Dormand-Prince 5(4) order method.
    // The Butcher Tableau is :
    constexpr long double c2 = 1.0 / 5.0;
    constexpr long double c3 = 3.0 / 10.0;
    constexpr long double c4 = 4.0 / 5.0;
    constexpr long double c5 = 8.0 / 9.0;
    constexpr long double a21 = 1.0 / 5.0;
    constexpr long double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    constexpr long double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    constexpr long double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    constexpr long double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    constexpr long double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

    constexpr long double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    constexpr long double e1 = 5179.0 / 57600.0, e3 = 7571.0 / 16695.0, e4 = 393.0 / 640.0, e5 = -92097.0 / 339200.0, e6 = 187.0 / 2100.0, e7 = 1.0 / 40.0;
}