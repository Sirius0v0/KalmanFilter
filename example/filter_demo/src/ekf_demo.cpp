#include <filter_demo/ekf_demo.h>
#include <cmath>

MyEKFilter::MyEKFilter() :EKFilter(9, 3, 4)
{
    this->Ts = 2;
    this->alpha = 1 / 10.0;
    this->SatAPos = { 1300e3, 400e3, 100e3 };
    this->SatBPos = { 1300e3, 400e3, -100e3 };
    this->acc_forecast.resize(this->get_u_dim(), 1);
    this->a_max.resize(this->get_u_dim(), 1);
    this->a_max.assign({ 2.0, 0.4, 0.8 });
}

MyEKFilter::MyEKFilter(int state_dim_, int u_dim_, int measure_dim_)
    :EKFilter(state_dim_, u_dim_, measure_dim_)
{
    this->Ts = 2;
    this->alpha = 1 / 10.0;
    this->SatAPos = { 1300e3, 400e3, 100e3 };
    this->SatBPos = { 1300e3, 400e3, -100e3 };
    this->acc_forecast.resize(this->get_u_dim(), 1);
    this->a_max.resize(this->get_u_dim(), 1);
    this->a_max.assign({ 2.0, 0.4, 0.8 });
}

void MyEKFilter::setA()
{
    kalmans::Matrix<double> A_onedim(3, 3);
    A_onedim.assign({
        1, Ts, (alpha * Ts - 1 + exp(-alpha * Ts)) / alpha / alpha,
        0, 1, (1 - exp(-alpha * Ts)) / alpha,
        0, 0, -exp(-alpha * Ts)
        });
    auto tmp = kalmans::Matrix<double>::BlkDiag(std::move(A_onedim), 3);
    this->A.swap(tmp);
}

void MyEKFilter::setB()
{
    kalmans::Matrix<double> B_onedim(3, 1);
    B_onedim.assign({
        (-Ts + alpha * Ts * Ts / 2 + (1 - exp(-alpha * Ts)) / alpha),
        Ts - (1 - exp(-alpha * Ts)) / alpha,
        1 - exp(-alpha * Ts)
        });
    auto tmp = kalmans::Matrix<double>::BlkDiag(std::move(B_onedim), 3);
    this->B.swap(tmp);
}

void MyEKFilter::setH()
{
    auto H_A = compute_H(this->SatAPos);
    auto H_B = compute_H(this->SatBPos);
    auto H_tmp = kalmans::Matrix<double>::ver_mat_cat({ std::move(H_A), std::move(H_B) });
    this->H.swap(H_tmp);
}

void MyEKFilter::setQ()
{
    double q11 = (2.0 * pow(alpha, 3) * pow(Ts, 3) - 6.0 * pow(alpha, 2) * pow(Ts, 2) + 6.0 * alpha * Ts + 3.0 - 12.0 * alpha * Ts * exp(-alpha * Ts) - 3.0 * exp(-2 * alpha * Ts)) / (6.0 * pow(alpha, 5));
    double q12 = (alpha * alpha * Ts * Ts - 2.0 * alpha * Ts + 1.0 - 2.0 * (1 - alpha * Ts) * exp(-alpha * Ts) + exp(-2 * alpha * Ts)) / (2.0 * pow(alpha, 4));
    double q13 = (1.0 - 2 * alpha * Ts * exp(-alpha * Ts) - exp(-2 * alpha * Ts)) / (2.0 * pow(alpha, 3));
    double q21 = q12;
    double q22 = (2 * alpha * Ts - 3 + 4 * exp(-alpha * Ts) - exp(-2 * alpha * Ts)) / (2.0 * pow(alpha, 3));
    double q23 = (1 - 2 * exp(-alpha * Ts) + exp(-2 * alpha * Ts)) / (2.0 * alpha * alpha);
    double q31 = q13;
    double q32 = q23;
    double q33 = (1 - exp(-2 * alpha * Ts)) / (2 * alpha);

    std::array<double, 3> sigma2{0.0};
    for (int i = 0; i < 3; ++i)
    {
        const double PI = 3.141592653589793;
        if (acc_forecast(i) >= 0)
            sigma2[i] = (4 - PI) / PI * pow((a_max(i) - acc_forecast(i)), 2);
        else
            sigma2[i] = (4 - PI) / PI * pow((a_max(i) + acc_forecast(i)), 2);
    }

    kalmans::Matrix tmp_qk{3, 3, { q11,q12,q13,q21,q22,q23,q31,q32,q33 }};
    tmp_qk *= (2 * alpha);
    auto tmp = kalmans::Matrix<double>::BlkDiag({ tmp_qk * sigma2[0], tmp_qk * sigma2[1],tmp_qk * sigma2[2] });
    this->Q.swap(tmp);
}

void MyEKFilter::setR()
{
    double err_std2 = pow(0.003 / 57.29577951308232, 2);
    this->R = kalmans::Matrix<double>::Diag({
        err_std2, err_std2, err_std2, err_std2
        });
}

kalmans::Matrix<double> MyEKFilter::compute_H(const std::array<double, 3>& satpos)
{
    double e1 = -((satpos[1] - x(3)) * (2 * satpos[0] - 2 * x(0)))
        / (2 * pow(
            (pow((satpos[0] - x(0)), 2)
                + pow((satpos[2] - x(6)), 2))
            , (3 / 2.0))
            * (pow((satpos[1] - x(3)), 2)
                / (pow((satpos[0] - x(0)), 2)
                    + pow((satpos[2] - x(6)), 2)) + 1));
    double e2 = pow(
        (pow((satpos[0] - x(0)), 2)
            + pow((satpos[2] - x(6)), 2))
        , (1 / 2.0))
        / (satpos[1] * satpos[1] - 2 * satpos[1] * x(3)
            + x(3) * x(3) +
            pow((satpos[0] - x(0)), 2)
            + pow((satpos[2] - x(6)), 2));
    double e3 = -((satpos[1] - x(3)) * (2 * satpos[2] - 2 * x(6)))
        / (2 * pow(
            (pow((satpos[0] - x(0)), 2)
                + pow((satpos[2] - x(6)), 2))
            , (3 / 2.0))
            * (pow((satpos[1] - x(3)), 2)
                / (pow((satpos[0] - x(0)), 2)
                    + pow((satpos[2] - x(6)), 2)) + 1));
    double e4 = (satpos[2] - x(6)) /
        (pow((satpos[0] - x(0)), 2) *
            (pow((satpos[2] - x(6)), 2) /
                pow((satpos[0] - x(0)), 2) + 1));
    double e5 = -1.0 / (
        (satpos[0] - x(0)) *
        (pow((satpos[2] - x(6)), 2)
            / pow((satpos[0] - x(0)), 2) + 1)
        );
    return kalmans::Matrix<double>(2, 9, {
        e1, 0, 0, e2, 0, 0, e3, 0, 0,
        e4, 0, 0, 0,  0, 0, e5, 0, 0
        });
}

void MyEKFilter::make_predict()
{
    this->acc_forecast.assign({ x(2),x(5),x(8) });
    this->x = this->A * this->x + this->B * this->acc_forecast;
}

void MyEKFilter::make_measure()
{
    this->z_hat(0) = atan((x(3) - SatAPos[1]) / sqrt(
        pow((x(0) - SatAPos[0]), 2) + pow((x(6) - SatAPos[2]), 2)));
    this->z_hat(1) = atan((x(6) - SatAPos[2]) / (x(0) - SatAPos[0]));
    this->z_hat(2) = atan((x(3) - SatBPos[1]) / sqrt(
        pow((x(0) - SatBPos[0]), 2) + pow((x(6) - SatBPos[2]), 2)));
    this->z_hat(3) = atan((x(6) - SatBPos[2]) / (x(0) - SatBPos[0]));
}

