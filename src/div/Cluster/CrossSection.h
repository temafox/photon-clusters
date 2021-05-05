#ifndef CLUSTER_DIV_CROSSSECTION_H
#define CLUSTER_DIV_CROSSSECTION_H

#include "Strip.h"

namespace cluster_div {

class CrossSection {
public:
    using strip_type = Strip;
    using id_type = strip_type::id_type;
    using float_type = double;

public:
    CrossSection() = delete;
    CrossSection(const id_type packedIDs[2], const float_type position[3]);
    CrossSection(id_type packedID0, id_type packedID1,
                 float_type x, float_type y, float_type z);

    inline const id_type *getPackedIDs() const noexcept { return packedIDs_; }
    inline id_type getLesserPackedID() const noexcept { return packedIDs_[0]; }
    inline id_type getGreaterPackedID() const noexcept { return packedIDs_[1]; }

    inline const float_type *getPosition() const noexcept { return position_; }
    inline float_type getX() const noexcept { return x_; }
    inline float_type getY() const noexcept { return y_; }
    inline float_type getZ() const noexcept { return z_; }

private:
    // Condition: [0] < [1]
    id_type packedIDs_[2] = {};

    float_type position_[3] = {};
    float_type &x_ = position_[0];
    float_type &y_ = position_[1];
    float_type &z_ = position_[2];
};

}

#endif //CLUSTER_DIV_CROSSSECTION_H
