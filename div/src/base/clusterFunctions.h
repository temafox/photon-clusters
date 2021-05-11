#ifndef CLUSTER_DIV_CLUSTERFUNCTIONS_H
#define CLUSTER_DIV_CLUSTERFUNCTIONS_H

#include "geraData.h"
#include "simData.h"
#include "functions.h"

namespace cluster_div {

ClusterMap *findClusterCenters(const gera_nm::strip_data &strips, const gera_nm::cross_data &cross_pos);

// Implementation

ClusterMap *findClusterCenters(const gera_nm::strip_data &strips, const gera_nm::cross_data &cross_pos) {
    auto *clusters = new ClusterMap();
    size_t numberOfCrosses = cross_pos.id1.size();

    class Running {
    public:
        double x, y, z;
        size_t size;

        Running() : x(0.), y(0.), z(0.), size(0) {}
    };
    std::map<Cluster_id_t, Running> running;

    for (size_t cross_i = 0; cross_i < numberOfCrosses; ++cross_i) {
        // Find the indices of the crossing strips
        size_t strip1 = findIndex(strips.packedID, cross_pos.id1[cross_i]);
        size_t strip2 = findIndex(strips.packedID, cross_pos.id2[cross_i]);

        // Ensure that either cluster_id != -1 (means no cluster)
        int cluster_id1 = strips.cluster_id[strip1];
        int cluster_id2 = strips.cluster_id[strip2];
        if ((cluster_id1 == -1) || (cluster_id2 == -1))
            continue;
        Cluster_id_t cluster_id_paired = Cluster_id_t(cluster_id1, cluster_id2);

        // Ensure that both clusters are on the same layer
        int layer1 = strips.layer[strip1];
        int layer2 = strips.layer[strip2];
        if (layer1 != layer2)
            continue;

        // Ensure that the current cluster pair is listed
        auto cl_it = clusters->find(cluster_id_paired);
        auto run_it = running.find(cluster_id_paired);
        if (cl_it == clusters->cend()) {
            // If not found, add the pair to the list
            cl_it = (clusters->insert(std::pair<Cluster_id_t, Cluster>(cluster_id_paired, Cluster()))).first;
            cl_it->second.layer = layer1;
            run_it = (running.insert(std::pair<Cluster_id_t, Running>(cluster_id_paired, Running()))).first;
        }

        run_it->second.x += cross_pos.x[cross_i];
        run_it->second.y += cross_pos.y[cross_i];
        run_it->second.z += cross_pos.z[cross_i];
        run_it->second.size += 1;
        cl_it->second.hitNumber += 1;
    }

    auto cl_it = clusters->begin();
    auto run_it = running.begin();
    for (; cl_it != clusters->end() && run_it != running.cend(); ++cl_it, ++run_it) {
        run_it->second.x /= run_it->second.size;
        run_it->second.y /= run_it->second.size;
        run_it->second.z /= run_it->second.size;

        // tg(phi) = y/x
        // tg(theta) = r/z, r^2 = x^2 + y^2
        cl_it->second.cphi = std::atan2(run_it->second.y, run_it->second.x);
        cl_it->second.ctheta = std::atan2(std::sqrt(std::pow(run_it->second.x, 2) + std::pow(run_it->second.y, 2)),
                                          run_it->second.z);
    }

    return clusters;
}

}

#endif // CLUSTER_DIV_CLUSTERFUNCTIONS_H
