#include <algorithm>
#include <stdexcept>

#include "LayerCluster.h"

namespace cluster_div {

const gera_nm::strip_data *LayerCluster::eventStrips_ = nullptr;
const gera_nm::cross_data *LayerCluster::eventCrossSections_ = nullptr;

LayerCluster::LayerCluster(std::size_t clusterID, const gera_nm::strip_data &stripData, const gera_nm::cross_data &crossData) {
    loadEvent(stripData, crossData);
    initLayerCluster_(clusterID);
}

LayerCluster::LayerCluster(std::size_t clusterID) {
    if (!eventStrips_ || !eventCrossSections_)
        throw std::logic_error{"Error in LayerCluster ctor: event uninitialized"};
    initLayerCluster_(clusterID);
}

void LayerCluster::loadEvent(const gera_nm::strip_data &stripData, const gera_nm::cross_data &crossData) {
    eventStrips_ = &stripData;
    eventCrossSections_ = &crossData;
}

void LayerCluster::initLayerCluster_(std::size_t clusterID) {
    auto stripIterator = std::find(eventStrips_->cluster_id.cbegin(), eventStrips_->cluster_id.cend(), clusterID);

    while (stripIterator != eventStrips_->cluster_id.cend()) {
        auto stripIndex = std::distance(eventStrips_->cluster_id.cbegin(), stripIterator);

        strips_.emplace_back(std::make_shared<strip_type>());
    }
}

}
