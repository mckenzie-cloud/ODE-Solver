#ifndef DENSE_OUT_H
#define DENSE_OUT_H

#include <vector>
#include <functional>

#include "../header/solver.h"

class DenseOut : Solver{
    public:
        DenseOut(std::function<void(long double, std::vector<long double> &, std::vector<long double> &)> fName, std::vector<std::vector<long double>> &solData, std::vector<long double> &tData, size_t dim);
        std::vector<long double> dense_eval_sol(long double t_val);
        void printDebug();
    private:
        std::vector<std::vector<long double>> solData;            // accumulated solution steps
        std::vector<long double> tData, m_K8, m_K9;               // accumulated time steps
        void process_two_extra_dopri_steps();
        std::vector<long double> interpolate(long double theta);
};

#endif