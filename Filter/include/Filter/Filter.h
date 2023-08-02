#pragma once

#include <Matrix/Matrix.hpp>

namespace kalmans
{
    class FilterBase
    {
    public:
        FilterBase();
        FilterBase(int state_dim_, int u_dim_, int measure_dim_);
        FilterBase(const FilterBase& rhs) = delete;
        FilterBase& operator=(const FilterBase& rhs) = delete;
        FilterBase(FilterBase&& rhs) = delete;
        virtual ~FilterBase() = default;

        void init(Matrix<double>& x, Matrix<double>& P);
        void step(Matrix<double>& u, Matrix<double>& z);
        const Matrix<double>& get_state() const;

        int get_state_dim() const;
        int get_u_dim() const;
        int get_measure_dim() const;

        void set_dim(int state_dim_, int u_dim_, int measure_dim_);

    protected:
        int state_dim;
        int u_dim;
        int measure_dim;

        Matrix<double> x;
        Matrix<double> u;
        Matrix<double> z;

        Matrix<double> A;
        Matrix<double> B;
        Matrix<double> H;
        Matrix<double> Q;
        Matrix<double> R;

        Matrix<double> z_hat;
        Matrix<double> P;
        Matrix<double> K;

        virtual void make_predict() = 0;
        virtual void time_update() = 0;
        virtual void make_measure() = 0;
        virtual void measure_update() = 0;

        void time_update_step(Matrix<double>& u);
        void measure_update_step(Matrix<double>& z);

        virtual void setA() {}
        virtual void setB() {}
        virtual void setH() {}
        virtual void setQ() {}
        virtual void setR() {}

    };
}
