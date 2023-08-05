#include <Filter/Filter.h>
#include <Matrix/Matrix.hpp>

namespace kalmans
{
    FilterBase::FilterBase()
    {
        this->state_dim = 0;
        this->u_dim = 0;
        this->measure_dim = 0;
        this->x.resize(state_dim, 1);
        this->u.resize(u_dim, 1);
        this->z.resize(measure_dim, 1);
        this->z_hat.resize(measure_dim, 1);
    }

    FilterBase::FilterBase(int state_dim_, int u_dim_, int measure_dim_)
    {
        this->state_dim = state_dim_;
        this->u_dim = u_dim_;
        this->measure_dim = measure_dim_;
        this->x.resize(state_dim, 1);
        this->u.resize(u_dim, 1);
        this->z.resize(measure_dim, 1);
        this->z_hat.resize(measure_dim, 1);
    }

    void FilterBase::set_dim(int state_dim_, int u_dim_, int measure_dim_)
    {
        this->state_dim = state_dim_;
        this->u_dim = u_dim_;
        this->measure_dim = measure_dim_;
    }

    int FilterBase::get_measure_dim() const
    {
        return this->measure_dim;
    }

    int FilterBase::get_state_dim() const
    {
        return this->state_dim;
    }

    int FilterBase::get_u_dim() const
    {
        return this->u_dim;
    }

    const Matrix<double>& FilterBase::get_state() const
    {
        return this->x;
    }

    void FilterBase::init(Matrix<double>& x, Matrix<double>& P)
    {
        this->x.swap(x);
        this->P.swap(P);
    }

    void FilterBase::step(Matrix<double>& u, Matrix<double>& z)
    {
        time_update_step(u);
        measure_update_step(z);
    }

    void FilterBase::time_update_step(Matrix<double>& u)
    {
        this->u.swap(u);

        this->setA();
        this->setB();
        this->make_predict();
        this->setQ();

        this->time_update();
    }

    void FilterBase::measure_update_step(Matrix<double>& z)
    {
        this->z.swap(z);

        this->setH();
        this->make_measure();
        this->setR();

        this->measure_update();
    }

}
