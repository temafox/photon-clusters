#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TF1.h"

#include <iostream>
#include <cmath>

int main( int argc, char *argv[] ) {
	/// Greeting
	std::cout << "Photon clusters - make data" << std::endl;
	if( argc < 2 ) {
		std::cerr << "Usage: " << argv[0] << " N" << std::endl;
		return -1;
	}

	/// Business part
	TRandom3 *rndm = new TRandom3( 0 );
	std::cout << "\nSeed: " << rndm->GetSeed() << std::endl;

	TF1 *gausHeight = new TF1( "gaus", "gaus(0)", 150, 200 );
	gausHeight->SetParameters( 1, 170, 20 );
	TF1 *gausMass = new TF1( "gaus", "gaus(0)", 0.5, 1.5 );
	gausMass->SetParameters( 1, 1, 0.1 );

	Float_t age, mass, height;
	Int_t type; // 1 = signal, -1 = background

	TFile file( "data/data.root", "RECREATE", "File with signal and background" );
	TTree *tr = new TTree( "tr", "signal = 1 & background = -1" );
	tr->Branch( "age", &age );
	tr->Branch( "mass", &mass );
	tr->Branch( "height", &height );
	tr->Branch( "type", &type );

	Int_t N = atoi( argv[argc-1] );
	for( Int_t i = 0; i < N; ++i ) {
		age = std::floor( rndm->Rndm() * 10 + 20 );
		
		// Adding a signal event
		if( rndm->Rndm() < 0.5 ) {
			Double_t p = {175};
			gausHeight->SetParameters( 1, p, 5 );

			height = gausHeight->GetRandom() * ( 1 - 0.002 * pow( age - 27, 2 ) );
			mass = ( height - 90 ) * gausMass->GetRandom();
			type = 1;
		}
		// Adding a background event
		else {
			Double_t p = {165};
			gausHeight->SetParameters( 1, p, 5 );

			height = gausHeight->GetRandom() * ( 1 - 0.002 * pow( age - 25, 2 ) );
			mass = ( height - 100 ) * gausMass->GetRandom();
			type = -1;
		}

		tr->Fill();
	}

	tr->Write();
	file.Close();

	return 0;
}
