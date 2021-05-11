#include "LXeEvent.h"

namespace cluster_div {

LXeEvent::LXeEvent(const gera_nm::tree_data &event)
:   nsim_{static_cast<uint_type>(event.nsim)},
    finalstate_id_{static_cast<uint_type>(event.finalstate_id)}
{
    for (auto i = 0U; i < nsim_; i++) {
        particles_.emplace_back(static_cast<ParticleType>(event.simtype[i]),
                                static_cast<ParticleType>(event.simorig[i]),
                                event.simmom[i], event.simphi[i], event.simtheta[i],
                                event.simvtx[i], event.simvty[i], event.simvtz[i]);
    }
}

}
