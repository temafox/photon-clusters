#include "sim_data.h"

namespace cluster_div {

	/// Create a tree of events where a photon comes from a pion decay
	void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" ) {
		TFile *inFile = new TFile( inFileName.c_str() );
		TTree *inTree = ( TTree * )inFile->Get( "tr_lxe" );

		Long64_t nEntries = inTree->GetEntries();

		sim_data event;
		gera_nm::strip_data *strips = new gera_nm::strip_data;
		gera_nm::cross_data *cross_pos = new gera_nm::cross_data;

		inTree->SetBranchAddress( "simtype", &(event.simtype) );
		inTree->SetBranchAddress( "simorig", &(event.simorig) );
		inTree->SetBranchAddress( "simmom", &(event.simmom) );
		inTree->SetBranchAddress( "simphi", &(event.simphi) );
		inTree->SetBranchAddress( "simtheta", &(event.simtheta) );
		inTree->SetBranchAddress( "simvtx", &(event.simvtx) );
		inTree->SetBranchAddress( "simvty", &(event.simvty) );
		inTree->SetBranchAddress( "simvtz", &(event.simvtz) );
		inTree->SetBranchAddress( "strips", &strips );
		inTree->SetBranchAddress( "cross_pos", &cross_pos );

		TFile *outFile = new TFile( outFileName.c_str(), "recreate" );
		TTree *outTree = inTree->CloneTree( 0 );
		
		for( Long64_t i = 0; i < nEntries; ++i ) {
			inTree->GetEntry();
			if( event.simtype == PHOTON && event.simorig == PION )
				outTree->Fill();
		}

		outFile->Print();
		delete inFile;
		delete outFile;
	}

}

void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" ) {
	cluster_div::selectPhotons( inFileName, outFileName );
}
