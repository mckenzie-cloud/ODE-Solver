#include "../header/dense_out.h"
#include "../header/solver.h"

DenseOut::DenseOut(std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName, 
                   std::vector<std::vector<long double>> &solData, std::vector<long double> &tData, size_t dim) : Solver(fName, dim), m_K8(dim), m_K9(dim)
{
    this->tData = tData;
    this->solData = solData;
}

std::vector<long double> DenseOut::dense_eval_sol(long double t_query)
{
    /**
     * Evaluating at theta = 0 gives yn solution (solData[i]) --while evaluating at theta = 1 gives yn_new solution (solData[i+1])
     */

    long double theta = 0.0L;
    for (size_t i = 0; i < tData.size() - 1; i++)
    {
        if (t_query >= tData[i] && t_query <= tData[i+1])
        {
            m_t  = tData[i];
            m_yn = solData[i];
            m_h  = tData[i+1] - tData[i];          // get the stepsize h

            process_doPri_steps();                 // calculate Kn steps at yn_i and h
            process_two_extra_dopri_steps();       // Process Two extra stages for 5th-order dense out interpolant.
            
            theta = (t_query - tData[i]) / m_h;    // theta ∈ [0, 1]
        }
    }
    return interpolate(theta);
}

void DenseOut::process_two_extra_dopri_steps()
{
    std::vector<long double> yn_temp(m_dim);

    for (size_t i = 0; i < m_dim; i++)
    {
        yn_temp[i] = m_yn[i] + m_h * (Dopri54::a81 * m_K1[i] + Dopri54::a83 * m_K3[i] +
                                      Dopri54::a84 * m_K4[i] + Dopri54::a85 * m_K5[i] + 
                                      Dopri54::a86 * m_K6[i] + Dopri54::a87 * m_K7[i]);
    }

    m_F(m_t + Dopri54::c8*m_h, yn_temp, m_K8);  //--------------------- 8th-stage -----------------------

    for (size_t i = 0; i < m_dim; i++)
    {
        yn_temp[i] = m_yn[i] + m_h * (Dopri54::a91 * m_K1[i] + Dopri54::a93 * m_K3[i] +
                                      Dopri54::a94 * m_K4[i] + Dopri54::a95 * m_K5[i] + 
                                      Dopri54::a96 * m_K6[i] + Dopri54::a97 * m_K7[i]);
    }

    m_F(m_t + Dopri54::c9*m_h, yn_temp, m_K9);  //--------------------- 9th-stage -----------------------

}

std::vector<long double> DenseOut::interpolate(long double t)  // 5th-order interpolant
{
    long double t_2 = t * t;
	long double t_3 = t_2 * t;
	long double t_4 = t_2 * t_2;

	long double b1 = (696.0L  *               t_4 - 2439.0L      * t_3 + 3104.0L   * t_2 - 1710.0L * t + 384.0L) /   384.0L;
	long double b3 = (500.0L  *  t * (24.0L * t_3 - 51.0L        * t_2 + 32.0L     * t   - 6.0L))                / (-1113.0L);
	long double b4 = (125.0L  *  t * (24.0L * t_3 - 51.0L        * t_2 + 32.0L     * t   - 6.0L))                / (-192.0L);
	long double b5 = (2187.0L *  t * (24.0L * t_3 - 51.0L        * t_2 + 32.0L     * t   - 6.0L))                /   6784.0L;
	long double b6 = (11.0L   *  t * (24.0L * t_3 - 51.0L        * t_2 + 32.0L     * t   - 6.0L))                / (-84.0L);
	long double b7 =            (t *         (t - 1.0L)          *      (32.0L     * t_2 - 31.0L   * t + 7.0L))  /   8.0L;
	long double b8 = (125.0L  *  t * std::pow(t - 1.0L, 2.0L))                                                   /   24.0L;
	long double b9 = (16.0L   *  t * std::pow(t - 1.0L, 2.0L)    *      (3.0L      * t   - 1.0L))                /   3.0L;

    std::vector<long double> solution(m_dim);
    for (size_t i = 0; i < m_dim; i++)
    {
        solution[i] = m_yn[i] + t * m_h * (b1 * m_K1[i] + b3 * m_K3[i] + b4 * m_K4[i] +
                                           b5 * m_K5[i] + b6 * m_K6[i] + b7 * m_K7[i] +
                                           b8 * m_K8[i] + b9 * m_K9[i]);
    }
    return solution;
}

void DenseOut::printDebug()
{
    for (size_t i = 0; i < m_dim; i++)
    {
        // display RK4 steps
        std::cout << m_K1[i] << '\n' <<
                     m_K2[i] << '\n' <<
                     m_K3[i] << '\n' <<
                     m_K4[i] << '\n' <<
                     m_K5[i] << '\n' <<
                     m_K6[i] << '\n' <<
                     m_K7[i] << '\n' <<
                     m_K8[i] << '\n' <<
                     m_K9[i] << '\n' <<
                     m_ynew[i] << '\n';
    }
    std::cout << '\n';
}