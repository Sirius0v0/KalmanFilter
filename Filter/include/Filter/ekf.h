#pragma once

#include <Filter/Filter.h>

namespace kalmans
{
    class EKFilter : public FilterBase
    {
    public:
        EKFilter() :FilterBase() {}
        EKFilter(int state_dim_, int u_dim_, int measure_dim_)
            :FilterBase(state_dim_, u_dim_, measure_dim_) {}
        ~EKFilter() override = default;

    protected:
        void time_update() override;
        void measure_update() override;

    private:
        Matrix<double> P_;
    };
}
