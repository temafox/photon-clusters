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
	std::map< Cluster_id_t, Cluster > const *findClusterCenters( gera_nm::strip_data const &strips, gera_nm::cross_data const &cross_pos );

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

		TH1F *angDist = new TH1F( "angDist", "Angular distance to closest cluster", 70, 0, TMath::Pi() );

		nEntries = 100;
		for( Long64_t i = 0; i < nEntries; ++i ) {
			inTree->GetEntry(i);
			
			std::map< Cluster_id_t, Cluster > const *clusters = findClusterCenters( *strips, *cross_pos );
			std::cout << "\nEntry #" << i << " processed" << std::endl;
			std::cout << "nsim = " << event.nsim << std::endl;

			for( int j = 0; j < event.nsim; ++j ) {
				std::cout << "simtype[" << j << "] = " << event.simtype[j]
					  << "\nsimorig[" << j << "] = " << event.simorig[j]
					  << "\nsimmom[" << j << "] = " << event.simmom[j]
					  << "\nsimphi[" << j << "] = " << event.simphi[j]
					  << "\nsimtheta[" << j << "] = " << event.simtheta[j]
					  << "\nsimvtx[" << j << "] = " << event.simvtx[j]
					  << "\nsimvty[" << j << "] = " << event.simvty[j]
					  << "\nsimvtz[" << j << "] = " << event.simvtz[j] << std::endl;

				if( event.simtype[j] == PHOTON && event.simorig[j] == PION ) {
					Photon photon(event, j);

					// For this photon find the closest cluster
					std::map< Cluster_id_t, Cluster >::const_iterator closestClusterIt = clusters->cbegin();
					std::cout << "clusters->size() = " << clusters->size() << std::endl;
					for( auto it = clusters->cbegin(); it != clusters->cend(); ++it ) {
						if( angularDistance(it->second, photon) < angularDistance(closestClusterIt->second, photon) )
							closestClusterIt = it;
					}
					if( closestClusterIt == clusters->cend() )
						continue;

					double angDistValue = angularDistance( closestClusterIt->second, photon );
					angDist->Fill( angDistValue );
				}
			}
		}
		
		std::cout << angDist->GetEntries() << std::endl;
		angDist->Draw();
	}

	std::map< Cluster_id_t, Cluster > const *findClusterCenters( gera_nm::strip_data const &strips, gera_nm::cross_data const &cross_pos ) {
		std::map< Cluster_id_t, Cluster > *clusters = new std::map< Cluster_id_t, Cluster >();
		size_t numberOfCrosses = cross_pos.id1.size();
		
		struct Running {
			double x, y, z;
			size_t size;

			Running(): x{0.}, y{0.}, z{0.}, size{0} {}
		};
		std::map< Cluster_id_t, Running > running;

		for( size_t cross_i = 0; cross_i < numberOfCrosses; ++cross_i ) {
			// Find the indices of the crossing strips
			size_t strip1 = findIndex( strips.packedID, cross_pos.id1[cross_i] );
			size_t strip2 = findIndex( strips.packedID, cross_pos.id2[cross_i] );

			// Ensure that either cluster_id != -1 (means no cluster)
			int cluster_id1 = strips.cluster_id[strip1];
			int cluster_id2 = strips.cluster_id[strip2];
			if( cluster_id1 == -1 || cluster_id2 == -1 )
				continue;
			Cluster_id_t cluster_id_paired = Cluster_id_t( cluster_id1, cluster_id2 );

			// Ensure that the current cluster pair is listed
			std::map< Cluster_id_t, Cluster >::iterator cl_it = clusters->find( cluster_id_paired );
			std::map< Cluster_id_t, Running >::iterator run_it = running.find( cluster_id_paired );
			if( cl_it == clusters->cend() ) { // not found
				cl_it = ( clusters->insert( std::pair< Cluster_id_t, Cluster >(cluster_id_paired, Cluster()) ) ).first;
				run_it = ( running.insert( std::pair< Cluster_id_t, Running >(cluster_id_paired, Running()) ) ).first;
			}

			run_it->second.x += cross_pos.x[cross_i];
			run_it->second.y += cross_pos.y[cross_i];
			run_it->second.z += cross_pos.z[cross_i];
			run_it->second.size += 1;
		}
		
		auto cl_it = clusters->begin();
		auto run_it = running.begin();
		for( ;
		     cl_it != clusters->end() && run_it != running.cend();
		     ++cl_it, ++run_it ) {
			run_it->second.x /= run_it->second.size;
			run_it->second.y /= run_it->second.size;
			run_it->second.z /= run_it->second.size;

			// tg(phi) = y/x
			// tg(theta) = r/z, r^2 = x^2 + y^2
			cl_it->second.cphi = std::atan2( run_it->second.y, run_it->second.x );
			cl_it->second.ctheta = std::atan2( std::sqrt( std::pow(run_it->second.x, 2) + std::pow(run_it->second.y, 2) ), run_it->second.z );
		}

		return clusters;
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
