#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>

#include "simData.h"

template <typename T>
size_t findIndex(const std::vector<T> &vector, const T &element);

template <typename T>
size_t findIndex(const T *array, size_t size, const T &element);

double angularDistance(const cluster_div::Angular &ang1, const cluster_div::Angular &ang2);
double phiDistance(double phi1, double phi2);

double energy(double mc2, double pc);

// Implementation

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

double angularDistance(const cluster_div::Angular &ang1, const cluster_div::Angular &ang2) {
    double phi1 = ang1.getPhi();
    double theta1 = ang1.getTheta();
    double phi2 = ang2.getPhi();
    double theta2 = ang2.getTheta();

    double angle = std::acos(std::cos(theta1) * std::cos(theta2)
        + std::sin(theta1) * std::sin(theta2) * std::cos(phi1 - phi2));
    return angle;
}

double phiDistance(double phi1, double phi2) {
    double diff = std::fabs(phi2 - phi1);
    return (diff > M_PI) ? (2. * M_PI - diff) : diff;
}

double energy(double mc2, double pc) {
    return std::sqrt(mc2 * mc2 + pc * pc);
}

#endif // FUNCTIONS_H
