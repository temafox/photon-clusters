#include <fstream>
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

void tmvaClassificationApplication( TString infilename = "data/data.root", TString outfilename = "data/tmva_app.root" ) {
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

	/// Create a histogram for events
	UInt_t nbin = 100;
	TH1F *histBDT(0);
	histBDT = new TH1F( "MVA_BDT", "Common MVA BDT response", nbin, -0.8, 0.8 );

	/// Open the input file
	TFile *inputFile(0);
	TString ifname = infilename;
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

	/// Extract the best cut value
	std::ifstream thresholdFile( "data/threshold.Float_t" );
	Float_t MVA_threshold;
	thresholdFile.read( ( char * )&MVA_threshold, sizeof( Float_t ) );
	thresholdFile.close();

	TFile *outputFile(0);
	TString ofname = outfilename;
	outputFile = new TFile( ofname, "RECREATE" );

	/// Set up the input tree
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;
	std::cout << "--- TMVAClassificationApplication: Using input file " << inputFile->GetName() << std::endl;
	std::cout << "--- Select signal sample" << std::endl;

	TTree *tr = ( TTree * )inputFile->Get( "tr" );
	tr->SetBranchAddress( "mass", &mass );
	tr->SetBranchAddress( "height", &height );
	tr->SetBranchAddress( "age", &age );
	Int_t type;
	tr->SetBranchAddress( "type", &type );

	/// Prepare the output trees
	TTree *tr_signal = new TTree( "tr_signal", "A tree of recognized signals" );
	TTree *tr_background = new TTree( "tr_background", "A tree of recognized background events" );
	
	tr_signal->Branch( "age", &age );
	tr_signal->Branch( "mass", &mass );
	tr_signal->Branch( "height", &height );
	tr_signal->Branch( "type", &type );

	tr_background->Branch( "age", &age );
	tr_background->Branch( "mass", &mass );
	tr_background->Branch( "height", &height );
	tr_background->Branch( "type", &type );

	/// Fill the histogram with the MVA function of events
	/// Discriminate events by the value of MVA
	std::cout << "MVA threshold: " << MVA_threshold << std::endl;

	Long64_t entryNumber = tr->GetEntries();
	std::cout << "--- Processing: " << entryNumber << " events" << std::endl;

	// Count false signal and background events
	Int_t falseBackgroundCount= 0, falseSignalCount= 0;

	TStopwatch stopwatch;
	stopwatch.Start();

	for( Long64_t i = 0; i < entryNumber; ++i ) {
		if( i % 1000 == 0 )
			std::cout << "--- ... Processing event: " << i << std::endl;
		tr->GetEntry( i );

		Float_t MVA_response = reader->EvaluateMVA( "BDT method" );
		histBDT->Fill( MVA_response );
		if( MVA_response < MVA_threshold ) {
			tr_background->Fill();
			if( type != -1 )
				++falseBackgroundCount;
		} else {
			tr_signal->Fill();
			if( type != 1 )
				++falseSignalCount;
		}
	}

	stopwatch.Stop();
	std::cout << "--- End of event loop: ";
	stopwatch.Print();

	/// Report on false recognitions
	std::cout << "Signal:"
		<< "\n- overall " << tr_signal->GetEntries()
		<< "\n- true    " << (tr_signal->GetEntries() - falseSignalCount)
		<< "\n- false   " << falseSignalCount<< std::endl;
	std::cout << "Background:"
		<< "\n- overall " << tr_background->GetEntries()
		<< "\n- true    " << (tr_background->GetEntries() - falseBackgroundCount)
		<< "\n- false   " << falseBackgroundCount<< std::endl;

	/// Write output, clean up
	histBDT->Write();
	tr_signal->Write();
	tr_background->Write();

	outputFile->Close();
	std::cout << "Created root file: " << outputFile->GetName() << " containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << std::endl;
}
