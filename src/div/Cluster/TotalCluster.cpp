#include "TotalCluster.h"

namespace cluster_div {

TotalCluster::TotalCluster(const LXeEvent &event) {
    ;
}

TotalCluster::layer_cluster_type &TotalCluster::atLayer(std::size_t layer, std::size_t clusterAtLayer) {
    return *layerClusters_[layer][clusterAtLayer];
}

const TotalCluster::layer_cluster_type &TotalCluster::atLayer(std::size_t layer, std::size_t clusterAtLayer) const {
    return *layerClusters_[layer][clusterAtLayer];
}

}