#include <algorithm>
#include <cmath>
#include "sim_data.h"

template < typename T >
size_t findIndex( std::vector< T > const &vec, T const &elem);
	
namespace cluster_div {

	/// Create a tree of events where a photon comes from a pion decay
	/// without saving into a file
	void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" );
	
	/// Create `Cluster' structures for each cluster
	std::vector< Cluster > const &findClusterCenters( gera_nm::strip_data const &strips, gera_nm::cross_data const &cross_pos );

	double angularDistance( Cluster const &, Photon const & );
	
	//////////////////////////////////////////////////////////////////
	//                        Implementations                       //
	//////////////////////////////////////////////////////////////////

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

		//TFile *outFile = new TFile( outFileName.c_str(), "recreate" );
		//TTree *outTree = inTree->CloneTree( 0 );
		
		TH1F *angDist = new TH1F( "angDist", "Angular distance to closest cluster", 70, 0, TMath::Pi() );

		nEntries = 1000;
		for( Long64_t i = 0; i < nEntries; ++i ) {
			inTree->GetEntry(i);
			
			std::vector< Cluster > clusters = findClusterCenters( *strips, *cross_pos );
			std::cout << "Entry #" << i << " processed" << std::endl;

			for( int j = 0; j < event.nsim; ++j ) {
				if( event.simtype[j] == PHOTON && event.simorig[j] == PION ) {
					Photon photon(event, j);

					// For this photon find the closest cluster
					size_t closestClusterID = -1;
					for( size_t clusterID = 0; clusterID < clusters.size(); ++clusterID ) {
						if( closestClusterID == -1 ) {
							closestClusterID = clusterID;
						}
						if( angularDistance(clusters[clusterID], photon) < angularDistance(clusters[closestClusterID], photon) )
							closestClusterID = clusterID;
					}
					if( closestClusterID == -1 )
						continue;

					std::cout << "clusters.size() = " << clusters.size() << "\nclosestClusterID = " << closestClusterID << std::endl;
					double angDistValue = angularDistance( clusters[closestClusterID], photon );
					angDist->Fill( angDistValue );
				}
			}
		}
		
		std::cout << angDist->GetEntries() << std::endl;
		angDist->Draw();

		//outFile->Print();
		delete inFile;
		//delete outFile;

		delete strips;
		delete cross_pos;
	}
	
	std::vector< Cluster > const &findClusterCenters( gera_nm::strip_data const &strips, gera_nm::cross_data const &cross_pos ) {
		std::vector< Cluster > *clusters = new std::vector< Cluster >();
		std::vector< double > runningX, runningY, runningZ;
		std::vector< size_t > sizes;
		size_t numberOfCrosses = cross_pos.id1.size();
		
		for( size_t cross_i = 0; cross_i < numberOfCrosses; ++cross_i ) {
			// Find the indices of the crossing strips
			size_t strip1 = findIndex( strips.packedID, cross_pos.id1[cross_i] );
			size_t strip2 = findIndex( strips.packedID, cross_pos.id2[cross_i] );

			// Ensure that both strips belong to the same cluster
			// and that cluster_id != -1 (means no cluster)
			int cluster_id = strips.cluster_id[strip1];
			if( cluster_id == -1 )
				continue;
			if( cluster_id != strips.cluster_id[strip2] )
				continue;

			// Find this cluster among the already listed ones
			int clusterIndex;
			std::vector< Cluster >::const_iterator clusterIt = std::find_if( clusters->cbegin(), clusters->cend(), [ cluster_id ]( Cluster const &c ) { return c.cluster_id == cluster_id; } );
			if( clusterIt == clusters->cend() ) {
				clusters->push_back( Cluster(cluster_id) );
				clusterIndex = clusters->size() - 1;
			} else {
				clusterIndex = clusterIt - clusters->cbegin();
			}

			// Add to the running sums
			runningX[clusterIndex] += cross_pos.x[cross_i];
			runningY[clusterIndex] += cross_pos.y[cross_i];
			runningZ[clusterIndex] += cross_pos.z[cross_i];
			++sizes[clusterIndex];
		}

		// Convert Cartesian into spherical
		for( size_t clusterIndex = 0; clusterIndex < clusters->size(); ++clusterIndex ) {
			runningX[clusterIndex] /= sizes[clusterIndex];
			runningY[clusterIndex] /= sizes[clusterIndex];
			runningZ[clusterIndex] /= sizes[clusterIndex];

			// tg(phi) = y/x
			// tg(theta) = r/z, r^2 = x^2 + y^2
			(*clusters)[clusterIndex].cphi = std::atan2( runningY[clusterIndex], runningX[clusterIndex] );
			(*clusters)[clusterIndex].ctheta = std::atan2( std::sqrt( std::pow(runningX[clusterIndex], 2) + std::pow(runningY[clusterIndex], 2) ), runningZ[clusterIndex] );

			std::cout << "cluster_id == " << (*clusters)[clusterIndex].cluster_id
			          << "\ncphi == " << (*clusters)[clusterIndex].cphi
				  << "\nctheta == " << (*clusters)[clusterIndex].ctheta
				  << std::endl;
		}

		return *clusters;
	}
	
	double angularDistance( Cluster const &cluster, Photon const &photon ) {
		double cphi = cluster.cphi;
		double ctheta = cluster.ctheta;
		double phphi = photon.simphi;
		double phtheta = photon.simtheta;

		return std::acos( std::cos(ctheta)*std::cos(phtheta) + std::sin(ctheta)*std::sin(phtheta)*std::cos(cphi-phphi) );
	}
}

template < typename T >
size_t findIndex( std::vector< T > const &vec, T const &elem) {
	typename std::vector< T >::const_iterator it = std::find( vec.cbegin(), vec.cend(), elem );
	if( it == vec.cend() ) // not found
		return -1;
	return it - vec.cbegin();
}

void selectPhotons( const std::string &inFileName = "/store25/semenov/strips_run039799.root", const std::string &outFileName = "/store25/bazhenov/piph.root" ) {
	cluster_div::selectPhotons( inFileName, outFileName );
}
