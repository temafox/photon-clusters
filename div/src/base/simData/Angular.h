#ifndef CLUSTER_DIV_ANGULAR_H
#define CLUSTER_DIV_ANGULAR_H

namespace cluster_div {

class Angular {
public:
    virtual ~Angular() = default;

    virtual double getPhi() const noexcept = 0;
    virtual double getTheta() const noexcept = 0;
};

}

#endif // CLUSTER_DIV_ANGULAR_H
