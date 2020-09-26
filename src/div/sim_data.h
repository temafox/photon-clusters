#ifndef SIM_DATA_H
#define SIM_DATA_H

#include <vector>
#include "TObject.h"

// Values from Monte Carlo simulation
#define PHOTON (22)
#define PION (111)

namespace gera_nm {
	
	struct strip_data {
		std::vector< int > packedID;
		std::vector< int > innerID;
		std::vector< int > layer;
		std::vector< int > direction;
		std::vector< int > cluster_id;
		std::vector< double > amp;
		std::vector< double > sigma;
	};
	
	struct cross_data {
		// id1 < id2, either is packedID
		std::vector< int > id1;
		std::vector< int > id2;
		std::vector< double > x;
		std::vector< double > y;
		std::vector< double > z;
	};

}

namespace cluster_div {

	struct sim_data {
		int simtype;
		int simorig;
		Float_t simmom;
		Float_t simphi;
		Float_t simtheta;
		Float_t simvtx;
		Float_t simvty;
		Float_t simvtz;
	};

}

#endif // SIM_DATA_H
