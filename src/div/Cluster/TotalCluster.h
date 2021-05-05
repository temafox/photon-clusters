#ifndef CLUSTER_DIV_TOTALCLUSTER_H
#define CLUSTER_DIV_TOTALCLUSTER_H

#include <memory>
#include <vector>

#include "LayerCluster.h"
#include "LXeTree/LXeEvent.h"

namespace cluster_div {

class TotalCluster {
public:
    using layer_cluster_type = LayerCluster;
    // By layer, then enumerating all the clusters on a layer
    using layer_cluster_container_type = std::vector<std::vector<std::shared_ptr<layer_cluster_type> > >;
    // Next layer, `npos` if last
    using layer_cluster_links_type = std::vector<std::vector<std::tuple<std::size_t, std::size_t> > >;

public:
    TotalCluster() = delete;
    explicit TotalCluster(const LXeEvent &event);

    layer_cluster_type &atLayer(std::size_t layer, std::size_t clusterAtLayer);
    const layer_cluster_type &atLayer(std::size_t layer, std::size_t clusterAtLayer) const;

private:
    layer_cluster_container_type layerClusters_;
    layer_cluster_links_type layerLinks_;
};

}

#endif //CLUSTER_DIV_TOTALCLUSTER_H
