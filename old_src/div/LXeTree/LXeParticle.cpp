#include "LXeParticle.h"

namespace cluster_div {

LXeParticle::LXeParticle(ParticleType simtype, ParticleType simorig,
                         const float_type momentum[3], const float_type vertex[3])
:   simtype_{simtype},
    simorig_{simorig}
{
    for (auto i = 0U; i < 3; i++) {
        momentum_[i] = momentum[i];
        vertex_[i] = vertex[i];
    }
}

LXeParticle::LXeParticle(ParticleType simtype, ParticleType simorig,
                         float_type simmom, float_type simphi, float_type simtheta,
                         float_type simvtx, float_type simvty, float_type simvtz)
:   simtype_{simtype},
    simorig_{simorig},
    momentum_{simmom, simphi, simtheta},
    vertex_{simvtx, simvty, simvtz}
{}

bool LXeParticle::isParticle(ParticleType type) const noexcept {
    return type == getSimtype();
}

bool LXeParticle::isFromParticle(ParticleType originType) const noexcept {
    return originType == getSimorig();
}

}