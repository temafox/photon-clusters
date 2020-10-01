#ifndef SIM_DATA_H
#define SIM_DATA_H

#include <vector>
#include "TObject.h"

// Values from Monte Carlo simulation
#define PHOTON (22)
#define PION (111)

#define MAX_SIM (100)

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

	struct tree_data {
		int nsim;
		int finalstate_id;
		int simtype[ MAX_SIM ];
		int simorig[ MAX_SIM ];
		Float_t simmom[ MAX_SIM ];
		Float_t simphi[ MAX_SIM ];
		Float_t simtheta[ MAX_SIM ];
		Float_t simvtx[ MAX_SIM ];
		Float_t simvty[ MAX_SIM ];
		Float_t simvtz[ MAX_SIM ];
	};

}

namespace cluster_div {
	
	struct Photon {
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
		Photon( gera_nm::tree_data const &, size_t index );
	};

/*
	struct Strip {
		
	};
	
	class Cross {
	private:
		int strip_id[2];
	public:
		Cross( int strip1, int strip2 ) {
			if( strip1 < strip2 ) {
				strip_id[0] = strip1;
				strip_id[1] = strip2;
			} else if( strip2 < strip1 ) {
				strip_id[0] = strip2;
				strip_id[1] = strip1;
			} else {
				throw std::exception(); // 
			}
		}
	};
*/
	
	struct Cluster {
		int cluster_id; // from gera_nm::strip_data

		double cphi;
		double ctheta;
		int numPhotons;
		
		Cluster( int cluster_id );
	};
	
}

///////////////////////////////////////////////////////////
//                    Implementations                    //
///////////////////////////////////////////////////////////

namespace cluster_div {

	int Photon::totalPhotonNumber = 0;

	Photon::Photon() {
		photonID = totalPhotonNumber;
		++totalPhotonNumber;
	}

	Photon::Photon( gera_nm::tree_data const &event, size_t index ) {
		this->simtype = event.simtype[index];
		this->simorig = event.simorig[index];
		this->simmom = event.simmom[index];
		this->simphi = event.simphi[index];
		this->simtheta = event.simtheta[index];
		this->simvtx = event.simvtx[index];
		this->simvty = event.simvty[index];
		this->simvtz = event.simvtz[index];
	}

	Cluster::Cluster( int cluster_id ):
		cluster_id{cluster_id},
		cphi{0.},
		ctheta{0.},
		numPhotons{0}
	{}
}
#endif // SIM_DATA_H
