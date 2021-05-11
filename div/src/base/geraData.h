#ifndef GERA_DATA_H
#define GERA_DATA_H

#include <TObject.h>

#include <vector>

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

#endif // GERA_DATA_H
