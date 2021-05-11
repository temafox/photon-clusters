#ifndef CLUSTER_DIV_CLUSTER_H
#define CLUSTER_DIV_CLUSTER_H

#include <vector>
#include <map>

#include "Angular.h"
#include "Photon.h"

namespace cluster_div {

class Cluster_id_t {
public:
    int cluster_id1;
    int cluster_id2;

    Cluster_id_t(int cluster_id1, int cluster_id2);

    bool operator<(const Cluster_id_t &right) const {
        return (cluster_id1 < right.cluster_id1) && (cluster_id2 < right.cluster_id2);
    }
};

class Cluster : public Angular {
public:
    int layer;
    double cphi;
    double ctheta;
    int hitNumber;
    std::vector<Photon> photons;

    int eventIndex;
    int pionParticleIndex;

    Cluster();

    double getPhi() const noexcept override { return cphi; }

    double getTheta() const noexcept override { return ctheta; }
};

typedef std::map<Cluster_id_t, Cluster> ClusterMap;

// Implementation

Cluster::Cluster() :
        layer(0),
        cphi(0.),
        ctheta(0.),
        hitNumber(0),
        photons(std::vector<Photon>()),
        eventIndex(0),
        pionParticleIndex(0) {}

Cluster_id_t::Cluster_id_t(int cluster_id1, int cluster_id2) :
        cluster_id1(std::min(cluster_id1, cluster_id2)),
        cluster_id2(std::max(cluster_id1, cluster_id2)) {}

}

#endif // CLUSTER_DIV_CLUSTER_H
