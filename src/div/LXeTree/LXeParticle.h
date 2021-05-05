#ifndef CLUSTER_DIV_LXEPARTICLE_H
#define CLUSTER_DIV_LXEPARTICLE_H

#include <cstdint>

#include <TMath.h>

namespace cluster_div {

enum class ParticleType : std::int32_t {
    NONE = 0,
    PHOTON = 22,
    PION_0 = 111
};

class LXeParticle {
public:
    using float_type = Float_t;

public:
    LXeParticle() = delete;
    LXeParticle(ParticleType simtype, ParticleType simorig,
                const float_type momentum[3], const float_type vertex[3]);
    LXeParticle(ParticleType simtype, ParticleType simorig,
                float_type simmom, float_type simphi, float_type simtheta,
                float_type simvtx, float_type simvty, float_type simvtz);

    bool isParticle(ParticleType type) const noexcept;
    bool isFromParticle(ParticleType originType) const noexcept;

    inline ParticleType getSimtype() const noexcept { return simtype_; }
    inline ParticleType getSimorig() const noexcept { return simorig_; }

    inline const float_type *getMomentum() const noexcept { return momentum_; }
    inline float_type getSimmom() const noexcept { return simmom_; }
    inline float_type getSimphi() const noexcept { return simphi_; }
    inline float_type getSimtheta() const noexcept { return simtheta_; }

    inline const float_type *getVertex() const noexcept { return vertex_; }
    inline float_type getSimvtx() const noexcept { return simvtx_; }
    inline float_type getSimvty() const noexcept { return simvty_; }
    inline float_type getSimvtz() const noexcept { return simvtz_; }

private:
    ParticleType simtype_;
    ParticleType simorig_;

    float_type momentum_[3] = {};
    float_type &simmom_ = momentum_[0];
    float_type &simphi_ = momentum_[1];
    float_type &simtheta_ = momentum_[2];

    float_type vertex_[3] = {};
    float_type &simvtx_ = vertex_[0];
    float_type &simvty_ = vertex_[1];
    float_type &simvtz_ = vertex_[2];
};

}

#endif //CLUSTER_DIV_LXEPARTICLE_H
