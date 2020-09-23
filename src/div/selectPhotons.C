const Int_t PHOTON = 22;
const Int_t PION = 111;

/// tr_lxe structures

typedef struct {
	Int_t nsim;
	std::vector< Int_t > simtype;
	std::vector< Int_t > simorig;
	std::vector< Float_t > simmom;
	std::vector< Float_t > simphi;
	std::vector< Float_t > simtheta;
	std::vector< Float_t > simvtx;
	std::vector< Float_t > simvty;
	std::vector< Float_t > simvtz;
} SimData;

typedef struct {
	std::vector< Int_t > packedID;
	std::vector< Int_t > innerID;
	std::vector< Int_t > layer;
	std::vector< Int_t > direction;
	std::vector< Int_t > cluster_id;
	std::vector< Float_t > amp;
	std::vector< Float_t > sigma;
} StripData;

typedef struct {
	// id1 < id2
	std::vector< Int_t > id1;
	std::vector< Int_t > id2;
	std::vector< Float_t > x;
	std::vector< Float_t > y;
	std::vector< Float_t > z;
} CrossData;

/// Create a tree of events where a photon comes from a pion decay
void selectPhotons( const std::string &fileName = "/store25/semenov/strips_run039799.root") {
	
	TFile *inputFile = new TFile(fileName.c_str());
	TTree *inputTree = ( TTree * )inputFile->Get( "tr_lxe;80" );
	
	Int_t photonCount = 0;
	
	SimData *simData = new SimData();
	Int_t nsim, simtype, simorig;
	Float_t simmom, simphi, simtheta;
	Float_t simvtx, simvty, simvtz;
	
	inputTree->SetBranchAddress( "nsim", &nsim );
	inputTree->SetBranchAddress( "simtype", &simtype );
	inputTree->SetBranchAddress( "simorig", &simorig );
	inputTree->SetBranchAddress( "simmom", &simmom );
	inputTree->SetBranchAddress( "simphi", &simphi );
	inputTree->SetBranchAddress( "simtheta", &simtheta );
	inputTree->SetBranchAddress( "simvtx", &simvtx );
	inputTree->SetBranchAddress( "simvty", &simvty );
	inputTree->SetBranchAddress( "simvtz", &simvtz );
	
	inputTree->GetEntry( 0 );
	for( Int_t i = 0; i < nsim; ++i ) {
		inputTree->GetEntry( i );

		if( simtype == PHOTON && simorig == PION ) {
			simData->simtype.push_back(simtype);
			simData->simorig.push_back(simorig);
			simData->simmom.push_back(simmom);
			simData->simphi.push_back(simphi);
			simData->simtheta.push_back(simtheta);
			simData->simvtx.push_back(simvtx);
			simData->simvty.push_back(simvty);
			simData->simvtz.push_back(simvtz);
			++photonCount;
		}
	}
	simData->nsim = photonCount;
	
	TTree *outputTree = new TTree( "tr_ph", "A tree for photons born in pi0 decay" );
	outputTree->
}
