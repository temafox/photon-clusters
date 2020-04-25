#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "TMVA/Types.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int tmvaClassification( TString infilename = "data/data.root", TString outfilename = "data/tmva.root" ) {
	/// Show a greeting
	std::cout << "Photon clusters - TMVA classification" << std::endl;

	/// Load the library
	TMVA::Tools::Instance();

	/// Open the input file
	TFile *inputFile(0);
	TString ifname = infilename;
	if( !gSystem->AccessPathName( ifname ) ) {
		inputFile = TFile::Open( ifname );
	} else {
		std::cerr << "ERROR: in src/algorithm/example/tmvaClassification.cpp" << std::endl;
		std::cerr << "*** Input file " << ifname << " not found" << std::endl;
		exit(1);
	}
	if( !inputFile ) {
		std::cerr << "ERROR: in src/algorithm/example/tmvaClassification.cpp" << std::endl;
		std::cerr << "*** Could not open input file " << ifname << std::endl;
		exit(1);
	}

	/// Open the output file <- create/overwrite it
	TString ofname = outfilename;
	TFile *outputFile = TFile::Open( ofname, "RECREATE" );

	/// Start classification
	std::cout << "==> Start TMVAClassification" << std::endl;
	std::cout << "--- TMVAClassification: Using input file " << inputFile->GetName() << std::endl;
	std::cout << "--- TMVAClassification: Using output file " << outputFile->GetName() << std::endl;

	/// Extract trees from the input file
	TTree *inputTree = ( TTree * )inputFile->Get( "tr" );
	TTree *signalTree = inputTree->CopyTree( "type == 1" );
	TTree *backgroundTree = inputTree->CopyTree( "type == -1" );
	
	/// Initialize a TMVA data loader
	TMVA::DataLoader *dataloader = new TMVA::DataLoader( "dataset" );
	dataloader->AddVariable( "mass", 'F' );
	dataloader->AddVariable( "height", 'F' );
	dataloader->AddVariable( "age", 'F' );

	Double_t signalWeight = 1.0;
	Double_t backgroundWeight = 1.0;

	dataloader->AddSignalTree( signalTree, signalWeight );
	dataloader->AddBackgroundTree( backgroundTree, backgroundWeight );

	TCut myCutSignal = "";
	TCut myCutBackground = "";

	dataloader->PrepareTrainingAndTestTree( myCutSignal, myCutBackground,
		"nTrain_Signal=0"
		":nTrain_Background=0"
		":SplitMode=Random"
		":NormMode=NumEvents"
		":!V" );

	/// Apply a TMVA factory to process the trees
	TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
		"!V"
		":!Silent"
		":Color"
		":DrawProgressBar"
		":Transformations=I;D;P;G,D"
		":AnalysisType=Classification" );

	factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
		"!H"
		":!V"
		":NTrees=850"
		":MinNodeSize=2.5%"
		":MaxDepth=3"
		":BoostType=AdaBoost"
		":AdaBoostBeta=0.5"
		":UseBaggedBoost"
		":BaggedSampleFraction=0.5"
		":SeparationType=GiniIndex"
		":nCuts=20" );

	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();

	/// Save the best cut value
	Float_t MVA_threshold = dynamic_cast< TMVA::MethodBase * >( factory->GetMethod( "dataset", "BDT" ) )->GetSignalReferenceCut();
	std::cout << "The MVA threshold value: " << MVA_threshold << std::endl;
	
	std::ofstream thresholdFile( "data/threshold.Float_t" );
	thresholdFile.write( ( char * )&MVA_threshold, sizeof( Float_t ) );
	thresholdFile << "\n" << MVA_threshold;

	/// Clean up
	thresholdFile.close();
	outputFile->Close();

	std::cout << "--- TMVAClassification: Wrote output" << std::endl;
	std::cout << "==> Finished TMVAClassification" << std::endl;

	delete factory;
	delete dataloader;

	// If we are not running in batch mode
	/// Launch a GUI for the output root macro
	if( !gROOT->IsBatch() )
		TMVA::TMVAGui( ofname );

	return 0;
}
