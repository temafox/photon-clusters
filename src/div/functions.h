#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>
#include "sim_data.h" 

template <typename T>
size_t findIndex(const std::vector<T> &vector, const T &element);

template <typename T>
size_t findIndex(const T *array, size_t size, const T &element);

double angularDistance(const cluster_div::Cluster &cluster, const cluster_div::Photon &photon);
double phiDistance(double phi1, double phi2);

double energy(double mc2, double pc);

namespace cluster_div {
ClusterMap *findClusterCenters(const gera_nm::strip_data &strips, const gera_nm::cross_data &cross_pos);
}

/// Implementations

template <typename T>
size_t findIndex(const std::vector<T> &vector, const T &element) {
    auto begin = vector.cbegin();
    auto end = vector.cend();

    auto elementIterator = std::find(begin, end, element);
    if (elementIterator == end) // element was not found
        return -1;
    return elementIterator - begin;
}

template <typename T>
size_t findIndex(const T *array, size_t size, const T &element) {
    for (size_t i = 0; i < size; ++i) {
        if (array[i] == element)
            return i;
    }
    return -1;
}

double angularDistance(const cluster_div::Cluster &cluster, const cluster_div::Photon &photon) {
    double clusterPhi = cluster.cphi;
    double clusterTheta = cluster.ctheta;
    double photonPhi = photon.simphi;
    double photonTheta = photon.simtheta;

    double angle = std::acos(std::cos(clusterTheta) * std::cos(photonTheta)
        + std::sin(clusterTheta) * std::sin(photonTheta) * std::cos(clusterPhi - photonPhi));
    return angle;
}

double phiDistance(double phi1, double phi2) {
    double diff = std::fabs(phi2 - phi1);
    return (diff > M_PI) ? (2. * M_PI - diff) : diff;
}

double energy(double mc2, double pc) {
    return std::sqrt(mc2 * mc2 + pc * pc);
}

namespace cluster_div {

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

#endif // FUNCTIONS_H
