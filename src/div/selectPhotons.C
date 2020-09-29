#include <algorithm>
#include <cmath>
#include "sim_data.h"

namespace cluster_div {

	/// Create a tree of events where a photon comes from a pion decay
	/// Without saving into a file
	void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" ) {
		TFile *inFile = new TFile( inFileName.c_str() );
		TTree *inTree = ( TTree * )inFile->Get( "tr_lxe;80" );

		Long64_t nEntries = inTree->GetEntries();

		gera_nm::tree_data event;
		gera_nm::strip_data *strips = new gera_nm::strip_data;
		gera_nm::cross_data *cross_pos = new gera_nm::cross_data;

		inTree->SetBranchAddress( "nsim", &(event.nsim) );
		inTree->SetBranchAddress( "simtype", event.simtype );
		inTree->SetBranchAddress( "simorig", event.simorig );
		inTree->SetBranchAddress( "simmom", event.simmom );
		inTree->SetBranchAddress( "simphi", event.simphi );
		inTree->SetBranchAddress( "simtheta", event.simtheta );
		inTree->SetBranchAddress( "simvtx", event.simvtx );
		inTree->SetBranchAddress( "simvty", event.simvty );
		inTree->SetBranchAddress( "simvtz", event.simvtz );
		inTree->SetBranchAddress( "strips", &strips );
		inTree->SetBranchAddress( "cross_pos", &cross_pos );

		TFile *outFile = new TFile( outFileName.c_str(), "recreate" );
		TTree *outTree = inTree->CloneTree( 0 );
		
		for( Long64_t i = 0; i < nEntries; ++i ) {
			inTree->GetEntry(i);
			
			std::vector< Cluster > clusters = findClusterCenters( strips, cross_pos );

			for( int j = 0; j < event.nsim; ++j ) {
				if( event.simtype[j] == PHOTON && event.simorig[j] == PION ) {
					// For this photon find the closest center of a cluster
					
					for( size_t clusterID = 0; clusterID < clusters.size(); ++clusterID ) {
						
					}
				}
			}
		}

		outFile->Print();
		delete inFile;
		delete outFile;

		delete strips;
		delete cross_pos;
	}
	
	std::vector< Cluster > const &findClusterCenters( gera_nm::strip_data const &strips, gera_nm::cross_data const &cross_pos ) {
		std::vector< Cluster > clusters;
		std::vector< double > runningX, runningY, runningZ;
		std::vector< size_t > sizes;
		size_t numberOfCrosses = cross_pos.id1.size();
		
		for( size_t cross_i = 0; cross_i < numberOfCrosses; ++cross_i ) {
			// Find the indices of the crossing strips
			std::vector< Cluster >::const_iterator strip1_iterator, strip2_iterator;
			strip1_iterator = std::find( strips.packedID.cbegin(), strips.packedID.cend(), cross_pos.id1[cross_i] );
			strip2_iterator = std::find( strips.packedID.cbegin(), strips.packedID.cend(), cross_pos.id2[cross_i] );

			size_t strip1 = strip1_iterator - strips.packedID.cbegin();
			size_t strip2 = strip2_iterator - strips.packedID.cbegin();

			// Ensure that both strips belong to the same cluster
			// and that clusterID != -1 (means no cluster)
			int clusterID = strips.cluster_id[strip1];
			if( clusterID == (size_t)(-1) )
				continue;
			if( clusterID != strips.cluster_id[strip2] )
				continue;

			// Expand `clusters' if we have come across a new cluster
			if( clusters.size() < clusterID + 1 ) {
				clusters.resize( clusterID + 1 );
				runningX.resize( clusterID + 1, 0 );
				runningY.resize( clusterID + 1, 0 );
				runningZ.resize( clusterID + 1, 0 );
				sizes.resize( clusterID + 1, 0 );
			}

			// Add to the running sums
			runningX[clusterID] += cross_pos.x[cross_i];
			runningY[clusterID] += cross_pos.y[cross_i];
			runningZ[clusterID] += cross_pos.z[cross_i];
			++sizes[clusterID];
		}

		// Convert Cartesian into spherical
		for( size_t clusterID = 0; clusterID < clusters.size(); ++clusterID ) {
			runningX[clusterID] /= sizes[clusterID];
			runningY[clusterID] /= sizes[clusterID];
			runningZ[clusterID] /= sizes[clusterID];

			// tg(phi) = y/x
			// tg(theta) = r/z, r^2 = x^2 + y^2
			clusters[clusterID].cphi = std::atan2( runningY[clusterID], runningX[clusterID] );
			clusters[clusterID].ctheta = std::atan2( std::sqrt( std::pow(runningX[clusterID], 2) + std::pow(runningY[clusterID], 2) ), runningZ[clusterID] );
		}
	}
	
}

void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" ) {
	cluster_div::selectPhotons( inFileName, outFileName );
}
