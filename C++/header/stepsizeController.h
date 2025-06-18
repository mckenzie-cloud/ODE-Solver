
namespace StepSizeController
{
    enum Controllers
    {
        STANDARD,
        H211PI,
        PI42,
        H211B,
        H312PID,
        H312B,
        H0321,
        H321

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
         * 5/4    | 1/2    |-3/4    |-1/4    |-3/4    | H0321    | smooth
         * 1/3    | 1/18   |-5/18   |-5/6    |-1/6    | H321     | medium
         *-------------------------------------------------------------------------
         */
    };
}