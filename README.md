# Kalman Filter
It is a generic implementation of Kalman Filter, should work for any system.

## Example

For the following filtering equation:

$$ \tilde{x}_{k+1}=f\left(\hat{x}_k, u_k, 0\right)+\Gamma_k\omega_k,\quad\omega_k\sim N(0,Q_k) $$ 

$$ \tilde{z}_{k+1}=h\left(\tilde{x}_{k+1}\right)+v_{k+1},\quad v_k\sim N(0,R_k) $$ 

linearize the discrete model around the current estimate, and the Jacobian matrix is as follows:

$$ \left. A_k=\frac{\partial f(x_k)}{\partial x_k}\right|_{x_k=\hat{x}_{k|k}} $$ 

$$ \left. B_k=\frac{\partial f(x_k)}{\partial u_k}\right|_{x_k=\hat{x}_{k|k}} $$ 

$$ \left. H_{k+1}=\frac{\partial h(x_{k+1})}{\partial x_{k+1}} \right|_{x_{k+1}=\hat{x}_{k+1|k}} $$ 

Here we will implement an Extended Kalman Filter to estimate the position, velocity, and acceleration of a target. The example can be found in `example/filter_demo`. You can access the complete demonstration there.

You need to inherit the `EKFilter` class, set the relevant parameters, and then override the state transition matrix and other Jacobian matrices in the filtering model:

```cpp
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
```

Perform initialization of the relevant parameters in the constructor function:

```cpp
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
```

Define the state transition matrix and other matrices:

```cpp
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
```

Override the state equation (`make_predict` function with the result assigned to `this->x_`) and the observation equation (`make_measure` function with the result assigned to `this->z_hat`):

```cpp
void MyEKFilter::make_predict()
{
    static bool isT0 = true;
    if (isT0)
    {
        this->acc_forecast.assign({ x(2),x(5),x(8) });
        isT0 = false;
    }
    else
    {
        this->acc_forecast.assign({ x_(2),x_(5),x_(8) });
    }
    this->x_ = this->A * this->x + this->B * this->acc_forecast;
}

void MyEKFilter::make_measure()
{
    this->z_hat(0) = atan((x_(3) - SatAPos[1]) / sqrt(
        pow((x_(0) - SatAPos[0]), 2) + pow((x_(6) - SatAPos[2]), 2)));
    this->z_hat(1) = atan((x_(6) - SatAPos[2]) / (x_(0) - SatAPos[0]));
    this->z_hat(2) = atan((x_(3) - SatBPos[1]) / sqrt(
        pow((x_(0) - SatBPos[0]), 2) + pow((x_(6) - SatBPos[2]), 2)));
    this->z_hat(3) = atan((x_(6) - SatBPos[2]) / (x_(0) - SatBPos[0]));
}
```

Set the initial estimated values and initial covariance, initialize the filter, and perform filtering based on the observations:

```cpp
int main()
{
    constexpr int STATE_DIM = 9;
    constexpr int U_DIM = 3;
    constexpr int MEASURE_DIM = 4;
    MyEKFilter filter(STATE_DIM, U_DIM, MEASURE_DIM);

    auto P0 = kalmans::Matrix<double>::Diag({
        1, 0.1, 0.01,
        1, 0.1, 0.01,
        1, 0.1, 0.01
        });
    kalmans::Matrix<double> x(STATE_DIM, 1);
    x(0) = 1500e3;
    x(3) = 10e3;
    x(7) = -250;

    int i = 0;
    // 1. initialize the filter
    filter.init(x, P0);

    kalmans::Matrix<double> z(MEASURE_DIM, 1);
    // When u is not required as input, simply declare it as a placeholder parameter
    kalmans::Matrix<double> u;
    filter_value.assign(filter.get_state(), STATE_DIM * i);

    for (; i < measure_data.get_row(); ++i)
    {
		// measure data
        for (int j = 0; j < MEASURE_DIM; ++j)
            z(j) = measure_data(i, j);
		// 2. Single-step filtering
        filter.step(u, z);
        // 3. get the filtered state data
        std::cout << "[Time=" << 2 * (i + 1) << "]" << filter.get_state().transpose() << '\n';
    }

    return EXIT_SUCCESS;
}
```

Plot the results:

![EKF滤波结果](https://cdn.jsdelivr.net/gh/Sirius0v0/image_store/blog/20230805153247.png)
