#ifndef CLUSTER_DIV_STRIP_H
#define CLUSTER_DIV_STRIP_H

#include <set>

namespace cluster_div {

class StripID {
public:
    using id_type = std::size_t;
    static constexpr std::size_t ID_COUNT = 5;

public:
    StripID() = delete;
    explicit StripID(const id_type stripIDInitializer[ID_COUNT]);
    StripID(id_type packedID, id_type innerID, id_type layer, id_type direction, id_type clusterID);

    bool operator<(const StripID &right) const;

    inline id_type getPackedID() const noexcept { return packedID_; }
    inline id_type getInnerID() const noexcept { return innerID_; }
    inline id_type getLayer() const noexcept { return layer_; }
    inline id_type getDirection() const noexcept { return direction_; }
    inline id_type getClusterID() const noexcept { return clusterID_; }

private:
    id_type IDs_[ID_COUNT] = {};

    id_type &packedID_ = IDs_[0];
    id_type &innerID_ = IDs_[1];
    id_type &layer_ = IDs_[2];
    id_type &direction_ = IDs_[3];
    id_type &clusterID_ = IDs_[4];
};

class Strip {
public:
    using strip_id_type = StripID;
    using id_type = strip_id_type::id_type;
    using float_type = double;

public:
    Strip() = delete;
    explicit Strip(id_type packedID);
    Strip(id_type packedID, float_type amplitude, float_type sigma);

    inline id_type getPackedID() const noexcept { return packedID_; }
    inline float_type getAmplitude() const noexcept { return amplitude_; }
    inline float_type getSigma() const noexcept { return sigma_; }

    inline void setAmplitude(float_type amplitude) { amplitude_ = amplitude; }
    inline void setSigma(float_type sigma) { sigma_ = sigma; }

private:
    id_type packedID_;

    float_type amplitude_ = 0.;
    float_type sigma_ = 0.;
};

}

#endif //CLUSTER_DIV_STRIP_H
