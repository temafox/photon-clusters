void searchFused( TString infilename = "data/data.root", TString outfilename = "data/searchFused.root" ) {
	/// Show a greeting
	std::cout << "Photon clusters - search for fused clusters" << std::endl;
	
	/// Open the input file
	TString ifname = infilename;
	TFile *inputFile(0);
	if( !gSystem->AccessPathName( ifname ) ) {
		inputFile = TFile::Open( ifname );
	} else {
		std::cerr << "ERROR: in src/algorithm/search/search_fused.cpp" << std::endl;
		std::cerr << "*** Input file " << ifname << " not found" << std::endl;
		exit(1);
	}
	if( !inputFile ) {
		std::cerr << "ERROR: in src/algorithm/search/search_fused.cpp" << std::endl;
		std::cerr << "*** Could not open input file " << ifname << std::endl;
		exit(1);
	}

	/// Prepare the input tree
	TTree *inputTree = ( TTree * )inputFile->Get( "tr_lxe;80" );
	Int_t nsim, simtype, simorig;
	Float_t simmom, simtheta, simphi;
	Float_t simvtx, simvty, simvtz;
	Float_t clusterr, clustertheta, clusterphi;

	inputTree->SetBranchAddress( "nsim", &nsim );
	inputTree->SetBranchAddress( "simtype", &simtype );
	inputTree->SetBranchAddress( "simorig", &simorig );
	inputTree->SetBranchAddress( "simmom", &simmom );
	inputTree->SetBranchAddress( "simtheta", &simtheta );
	inputTree->SetBranchAddress( "simphi", &simphi );
	inputTree->SetBranchAddress( "simvtx", &simvtx );
	inputTree->SetBranchAddress( "simvty", &simvty );
	inputTree->SetBranchAddress( "simvtz", &simvtz );
	
	inputTree->GetEntry( 298 );
	std::cout << "nsim = " << nsim
	          << "\nsimtype = " << simtype
	          << "\nsimorig = " << simorig
	          << "\nsimmom = " << simmom
	          << "\nsimtheta = " << simtheta
	          << "\nsimphi = " << simphi
	          << "\nsimvtx = " << simvtx
	          << "\nsimvty = " << simvty
	          << "\nsimvtz = " << simvtz;

	/// Clean up
	inputFile->Close();
}
