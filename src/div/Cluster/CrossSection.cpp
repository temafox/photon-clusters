#include "CrossSection.h"

namespace cluster_div {

CrossSection::CrossSection(const id_type packedIDs[2], const float_type position[3]) {
    packedIDs_[0] = std::min(packedIDs[0], packedIDs[1]);
    packedIDs_[1] = std::max(packedIDs[0], packedIDs[1]);

    x_ = position[0];
    y_ = position[1];
    z_ = position[2];
}

CrossSection::CrossSection(id_type packedID0, id_type packedID1,
        float_type x, float_type y, float_type z)
:   packedIDs_{std::min(packedID0, packedID1), std::max(packedID0, packedID1)},
    position_{x, y, z}
{}

}