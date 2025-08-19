
namespace Dopri54
{
    // We will be using the Dormand-Prince 5(4) order method.
    // The Butcher Tableau is :
    constexpr long double c2  = 1.0L     / 5.0L;
    constexpr long double c3  = 3.0L     / 10.0L;
    constexpr long double c4  = 4.0L     / 5.0L;
    constexpr long double c5  = 8.0L     / 9.0L;
    constexpr long double a21 = 1.0L     / 5.0L;
    constexpr long double a31 = 3.0L     / 40.0L,   a32 = 9.0L      / 40.0L;
    constexpr long double a41 = 44.0L    / 45.0L,   a42 = -56.0L    / 15.0L,    a43 = 32.0L    / 9.0L;
    constexpr long double a51 = 19372.0L / 6561.0L, a52 = -25360.0L / 2187.0L,  a53 = 64448.0L / 6561.0L, a54 = -212.0L   / 729.0L;
    constexpr long double a61 = 9017.0L  / 3168.0L, a62 = -355.0L   / 33.0L,    a63 = 46732.0L / 5247.0L, a64 = 49.0L     / 176.0L,    a65 = -5103.0L / 18656.0L;
    constexpr long double a71 = 35.0L    / 384.0L,                              a73 = 500.0L   / 1113.0L, a74 = 125.0L    / 192.0L,    a75 = -2187.0L / 6784.0L, a76 = 11.0L / 84.0L;

    constexpr long double b1  = 35.0L    / 384.0L,   b3 = 500.0L    / 1113.0L,  b4  = 125.0L   / 192.0L,  b5  = -2187.0L  / 6784.0L,   b6  = 11.0L    / 84.0L;
    constexpr long double e1  = 5179.0L  / 57600.0L, e3 = 7571.0L   / 16695.0L, e4  = 393.0L   / 640.0L,  e5  = -92097.0L / 339200.0L, e6  = 187.0L   / 2100.0L, e7  = 1.0L / 40.0L;

    // parameters for two-extra stages for 5th-order dense output interpolant for Dopri5.
    constexpr long double c8 = 1.0L      / 5.0L;
    constexpr long double c9 = 1.0L      / 2.0L;
    
    constexpr long double a81 = 5207.0L  / 48000.0L, a83 = 92.0L    / 795.0L,   a84 = -79.0L   / 960.0L,  a85 = 53217.0L  / 848000.0L, a86 = -11.0L   / 300.0L, a87 = 4.0L  / 125.0L;
    constexpr long double a91 = 613.0L   / 6144.0L,  a93 = 125.0L   / 318.0L,   a94 = -125.0L  / 3072.0L, a95 = 8019.0L   / 108544.0L, a96 = -11.0L   / 192.0L, a97 = 1.0L  / 32.0L;
}