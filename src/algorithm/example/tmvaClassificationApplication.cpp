#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

void tmvaClassificationApplication() {
	/// Show a greeting
	std::cout << "Photon clusters - TMVA classification application" << std::endl;

	/// Load the library
	TMVA::Tools::Instance();

	/// Create a reader to work with BDT
	TMVA::Reader *reader = new TMVA::Reader( 
			"!Color"
			":!Silent" );
	Float_t mass, height, age;

	reader->AddVariable( "mass", &mass );
	reader->AddVariable( "height", &height );
	reader->AddVariable( "age", &age );

	TString dir = "dataset/weights/";
	TString prefix = "TMVAClassification";
	TString methodName = TString( "BDT method" );
	TString weightFile = dir + prefix + TString( "_BDT.weights.xml" );

	reader->BookMVA( methodName, weightFile );

	/// Create a histogram for training events
	UInt_t nbin = 100;
	TH1F *histBDT(0);
	histBDT = new TH1F( "MVA_BDT", "MVA_BDT", nbin, -0.8, 0.8 );

	/// Open the input file
	TFile *inputFile(0);
	TString ifname = "data/data.root";
	if( !gSystem->AccessPathName( ifname ) ) {
		inputFile = TFile::Open( ifname );
	} else {
		std::cerr << "ERROR: in src/algorithm/example/tmvaApplication.cpp" << std::endl;
		std::cerr << "*** Input file " << ifname << " not found" << std::endl;
		exit(1);
	}
	if( !inputFile ) {
		std::cerr << "ERROR: in src/algorithm/example/tmvaApplication.cpp" << std::endl;
		std::cerr << "*** Could not open input file " << ifname << std::endl;
		exit(1);
	}

	/// Open the output file <- create/overwrite it
	TFile *outputFile(0);
	TString ofname = "data/tmva_app.root";
	outputFile = new TFile( ofname, "RECREATE" );

	/// Set up the input tree
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;
	std::cout << "--- TMVAClassificationApplication: Using input file " << inputFile->GetName() << std::endl;
	std::cout << "--- Select signal sample" << std::endl;

	TTree *tr = ( TTree * )inputFile->Get( "tr" );
	tr->SetBranchAddress( "mass", &mass );
	tr->SetBranchAddress( "height", &height );
	tr->SetBranchAddress( "age", &age );

	/// Fill the histogram with MVA function of events
	Long64_t entryNumber = tr->GetEntries();
	std::cout << "--- Processing: " << entryNumber << " events" << std::endl;

	TStopwatch stopwatch;
	stopwatch.Start();

	for( Long64_t i = 0; i < entryNumber; ++i ) {
		if( i % 1000 == 0 )
			std::cout << "--- ... Processing event: " << i << std::endl;
		tr->GetEntry( i );
		histBDT->Fill( reader->EvaluateMVA( "BDT method" ) );
	}

	stopwatch.Stop();
	std::cout << "--- End of event loop: ";
	stopwatch.Print();

	/// Write output, clean up
	histBDT->Write();

	outputFile->Close();
	std::cout << "Created root file: " << outputFile->GetName() << " containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << std::endl;
}
