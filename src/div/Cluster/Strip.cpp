#include "Strip.h"

namespace cluster_div {

StripID::StripID(const id_type stripIDInitializer[ID_COUNT])
{
    for (auto i = 0U; i < ID_COUNT; i++) {
        IDs_[i] = stripIDInitializer[i];
    }
}

StripID::StripID(id_type packedID, id_type innerID, id_type layer, id_type direction, id_type clusterID)
:   IDs_{packedID, innerID, layer, direction, clusterID}
{}

bool StripID::operator<(const StripID &right) const {
    return packedID_ < right.packedID_;
}

Strip::Strip(id_type packedID)
:   packedID_{packedID}
{}

Strip::Strip(id_type packedID, float_type amplitude, float_type sigma)
:   packedID_{packedID},
    amplitude_{amplitude},
    sigma_{sigma}
{}

}