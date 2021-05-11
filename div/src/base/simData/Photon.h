#ifndef CLUSTER_DIV_PHOTON_H
#define CLUSTER_DIV_PHOTON_H

#include <TObject.h>

#include "Angular.h"

namespace cluster_div {

class Photon : public Angular {
public:
    static int totalPhotonNumber;
    int photonID;

    int simtype;
    int simorig;

    Float_t simmom;
    Float_t simphi;
    Float_t simtheta;

    Float_t simvtx;
    Float_t simvty;
    Float_t simvtz;

    Photon();

    Photon(const gera_nm::tree_data &, size_t index);

    double getPhi() const noexcept override { return simphi; }

    double getTheta() const noexcept override { return simtheta; }
};

// Implementation

int Photon::totalPhotonNumber = 0;

Photon::Photon() :
        photonID(totalPhotonNumber++),
        simtype(PHOTON),
        simorig(0),
        simmom(0),
        simphi(0),
        simtheta(0),
        simvtx(0),
        simvty(0),
        simvtz(0) {}

Photon::Photon(const gera_nm::tree_data &event, size_t index) :
        photonID(totalPhotonNumber++),
        simtype(event.simtype[index]),
        simorig(event.simorig[index]),
        simmom(event.simmom[index]),
        simphi(event.simphi[index]),
        simtheta(event.simtheta[index]),
        simvtx(event.simvtx[index]),
        simvty(event.simvty[index]),
        simvtz(event.simvtz[index]) {}

}

#endif // CLUSTER_DIV_PHOTON_H
