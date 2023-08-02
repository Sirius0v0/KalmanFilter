#include <Filter/ekf.h>

namespace kalmans
{

    void EKFilter::time_update()
    {
        this->P_ = this->A * this->P * this->A.transpose() + this->Q;
    }

    void EKFilter::measure_update()
    {
        this->K = this->P_ * this->H.transpose() *
            (this->H * this->P * this->H.transpose() + this->R).inverse();
        this->x = this->x + this->K * (this->z - this->z_hat);
        this->P = (Matrix<double>::Eye(this->get_state_dim()) - this->K * this->H) * this->P_;
    }

}
