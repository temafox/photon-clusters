#ifndef CLUSTER_DIV_LAYERCLUSTER_H
#define CLUSTER_DIV_LAYERCLUSTER_H

#include <memory>
#include <vector>

#include "sim_data.h"
#include "Strip.h"
#include "CrossSection.h"

namespace cluster_div {

class LayerCluster {
public:
    using strip_type = Strip;
    using strip_container_type = std::vector<std::shared_ptr<strip_type> >;

    using crosssection_type = CrossSection;
    using crosssection_container_type = std::vector<std::shared_ptr<crosssection_type> >;

public:
    LayerCluster() = default;

    LayerCluster(std::size_t clusterID, const gera_nm::strip_data &stripData, const gera_nm::cross_data &cross_data);
    explicit LayerCluster(std::size_t clusterID);
    static void loadEvent(const gera_nm::strip_data &stripData, const gera_nm::cross_data &crossData);

    inline const strip_container_type &getStrips() const noexcept { return strips_; }
    inline strip_container_type &getStrips() { return strips_; }

    inline const crosssection_container_type &getCrossSections() const noexcept { return crossSections_; }
    inline crosssection_container_type &getCrossSections() { return crossSections_; }

private:
    void initLayerCluster_(std::size_t clusterID);

    static const gera_nm::strip_data *eventStrips_;
    static const gera_nm::cross_data *eventCrossSections_;

    strip_container_type strips_ = {};
    crosssection_container_type crossSections_ = {};
};

}

#endif //CLUSTER_DIV_LAYERCLUSTER_H
