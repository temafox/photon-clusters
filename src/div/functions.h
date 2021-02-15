#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include "sim_data.h"

template <typename T>
size_t findIndex(const std::vector<T> &vector, const T &element);

double angularDistance(const cluster_div::Cluster &cluster, const cluster_div::Photon &photon);

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

double angularDistance(const cluster_div::Cluster &cluster, const cluster_div::Photon &photon) {
    double clusterPhi = cluster.cphi;
    double clusterTheta = cluster.ctheta;
    double photonPhi = photon.simphi;
    double photonTheta = photon.simtheta;

    double angle = std::acos(std::cos(clusterTheta) * std::cos(photonTheta)
        + std::sin(clusterTheta) * std::sin(photonTheta) * std::cos(clusterPhi - photonPhi));
    return angle;
}

#endif // FUNCTIONS_H