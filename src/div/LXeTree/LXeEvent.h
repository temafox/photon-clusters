#ifndef CLUSTER_DIV_LXEEVENT_H
#define CLUSTER_DIV_LXEEVENT_H

#include <vector>
#include <cstdint>

#include "LXeParticle.h"
#include "sim_data.h"

namespace cluster_div {

class LXeEvent {
public:
    using uint_type = std::uint32_t;

public:
    LXeEvent() = delete;
    explicit LXeEvent(const gera_nm::tree_data &event);

private:
    uint_type nsim_;
    uint_type finalstate_id_;

    std::vector<LXeParticle> particles_;
};

}

#endif //CLUSTER_DIV_LXEEVENT_H
