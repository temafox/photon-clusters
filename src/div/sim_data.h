#ifndef SIM_DATA_H
#define SIM_DATA_H

#include <TObject.h>

#include <vector>
#include <cmath>

// Values from Monte Carlo simulation
#define PHOTON (22)
#define PION (111)

#define LAYERS (7)

#define MAX_SIM (100)

namespace gera_nm {
	
struct strip_data {
    std::vector<int> packedID;
    std::vector<int> innerID;
    std::vector<int> layer;
    std::vector<int> direction;
    std::vector<int> cluster_id;
    std::vector<double> amp;
    std::vector<double> sigma;
};

struct cross_data {
    // id1 < id2, both are packedID
    std::vector<int> id1;
    std::vector<int> id2;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

struct tree_data {
    int nsim;
    int finalstate_id;
    int simtype[MAX_SIM];
    int simorig[MAX_SIM];
    Float_t simmom[MAX_SIM];
    Float_t simphi[MAX_SIM];
    Float_t simtheta[MAX_SIM];
    Float_t simvtx[MAX_SIM];
    Float_t simvty[MAX_SIM];
    Float_t simvtz[MAX_SIM];
};

}

namespace cluster_div {

class Photon {
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
};

class Cluster_id_t {
public:
    int cluster_id1;
    int cluster_id2;

    Cluster_id_t(int cluster_id1, int cluster_id2);
    bool operator<(const Cluster_id_t &right) const {
        return (cluster_id1 < right.cluster_id1) && (cluster_id2 < right.cluster_id2);
    }
};

class Cluster {
public:
    int layer;
    double cphi;
    double ctheta;
    std::vector<Photon> photons;

    Cluster();
};

typedef std::map<Cluster_id_t, Cluster> ClusterMap;

class MultipleHistos {
public:
    TH1F **histArray;
    size_t size;
    std::string nameStart;
    std::string titleStart;
    size_t bins;
    double minX;
    double maxX;

    MultipleHistos(size_t size, std::string nameStart, std::string titleStart, size_t bins, double minX, double maxX);
    ~MultipleHistos();

    TH1F &operator[](size_t layer) const;
private:
    void initMultipleHistos();
};

}

/// Implementations

namespace cluster_div {

int Photon::totalPhotonNumber = 0;

Photon::Photon():
    photonID(totalPhotonNumber++)
{}

Photon::Photon(const gera_nm::tree_data &event, size_t index):
    simtype(event.simtype[index]),
    simorig(event.simorig[index]),
    simmom(event.simmom[index]),
    simphi(event.simphi[index]),
    simtheta(event.simtheta[index]),
    simvtx(event.simvtx[index]),
    simvty(event.simvty[index]),
    simvtz(event.simvtz[index])
{}

Cluster::Cluster():
    cphi(0.),
    ctheta(0.),
    photons(std::vector<Photon>())
{}

Cluster_id_t::Cluster_id_t(int cluster_id1, int cluster_id2):
    cluster_id1(std::min(cluster_id1, cluster_id2)),
    cluster_id2(std::max(cluster_id1, cluster_id2))
{}

MultipleHistos::MultipleHistos(size_t size, std::string nameStart, std::string titleStart, size_t bins, double minX, double maxX):
    size(size),
    nameStart(std::move(nameStart)),
    titleStart(std::move(titleStart)),
    bins(bins),
    minX(minX),
    maxX(maxX)
{
    histArray = new TH1F*[size];
    initMultipleHistos();
}

MultipleHistos::~MultipleHistos() {
    for (int i = 0; i < size; ++i)
        delete histArray[i];
    delete[] histArray;
}

TH1F &MultipleHistos::operator[](size_t layer) const {
    return *(histArray[layer]);
}

void MultipleHistos::initMultipleHistos() {
    histArray = new TH1F*[size];

    for (int i = 0; i < size; ++i) {
        std::string name(nameStart), title(titleStart);
        name.append(std::to_string(i));
        title.append(std::to_string(i));

        histArray[i] = new TH1F(name.c_str(), title.c_str(), bins, minX, maxX);
    }
}

}

#endif // SIM_DATA_H
