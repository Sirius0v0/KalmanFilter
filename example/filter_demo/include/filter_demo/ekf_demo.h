#pragma once

#include <Filter/ekf.h>
#include <array>

using kalmans::EKFilter;

class MyEKFilter : public EKFilter
{
public:
    MyEKFilter();
    MyEKFilter(int state_dim_, int u_dim_, int measure_dim_);
    ~MyEKFilter() override = default;

protected:
    void setA() override;
    void setB() override;
    void setH() override;
    void setQ() override;
    void setR() override;
    void make_predict() override;
    void make_measure() override;

private:
    double Ts;
    double alpha;
    std::array<double, 3> SatAPos;
    std::array<double, 3> SatBPos;

    kalmans::Matrix<double> acc_forecast;
    kalmans::Matrix<double> a_max;

    kalmans::Matrix<double> compute_H(const std::array<double, 3>& satpos);
};