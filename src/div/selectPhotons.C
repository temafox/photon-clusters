#include <algorithm>
#include <cmath>
#include "sim_data.h"

/// Search for `elem` in `vec` and return its index
template < typename T >
size_t findIndex(
    std::vector< T > const &vec,
    T const &elem
);

/// Construct a 2D array `histArray` of TH1F
/// with `size` histograms numbered 0..(`size`-1)
void initLayeredHistos(
    TH1F **histArray,
    size_t size,
    std::string const &nameStart,
    std::string const &titleStart,
    size_t bins,
    double minX,
    double maxX
);

namespace cluster_div {

    /// Create a tree of events where a photon comes from a pion decay
    /// without saving into a file
    void selectPhotons(
        const std::string &inFileName = "/store25/semenov/strips_run039799.root",
        const std::string &outFileName = "/store25/bazhenov/piph.root"
    );
    
    typedef std::map< Cluster_id_t, Cluster > Cluster_map;
    /// Syntactic sugar for `it->second`
    Cluster &cluster( Cluster_map::iterator const &it ) {
        return it->second;
    }
    Cluster const &cluster( Cluster_map::const_iterator const &it ) {
        return it->second;
    }

    /// Create `Cluster` structures for each cluster
    Cluster_map const *findClusterCenters(
        gera_nm::strip_data const &strips,
            gera_nm::cross_data const &cross_pos
    );

    /// Calculate the angular distance between the direction
    /// of the center of a cluster and the direction of a photon
    double angularDistance(
        Cluster const &cluster,
        Photon const &photon
    );

}
    
////////////////////////////////////////////////////////////////////////////////
//                              Implementations                               //
////////////////////////////////////////////////////////////////////////////////

namespace cluster_div {

    void selectPhotons(
        const std::string &inFileName = "/store25/semenov/strips_run039799.root",
        const std::string &outFileName = "/store25/bazhenov/piph.root"
    ) {
        // Initialize inputs
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

        // Create histograms of angular distances
        // Layers are numbered 0..6
        TH1F *angDist = new TH1F(
            "angDist",
            "Angular distance to closest cluster",
            70, 0, TMath::Pi()
        );

        TH1F *angDistLayer[LAYERS];
        TH1F *phiDistLayer[LAYERS];
        TH1F *thetaDistLayer[LAYERS];
        initLayeredHistos(
            angDistLayer,
            LAYERS,
            "angDist",
            "Angular distance at layer ",
            70, 0, TMath::Pi()
        );
        initLayeredHistos(
            phiDistLayer,
            LAYERS,
            "phiDist",
            "Difference for phi at layer ",
            140, 0, 2. * TMath::Pi()
        );
        initLayeredHistos(
            thetaDistLayer,
            LAYERS,
            "thetaDist",
            "Difference for theta at layer ",
            70, 0, TMath::Pi()
        );

        // Process events
        nEntries = 100;
        std::cout << "nEntries = " << nEntries << std::endl;
        for( Long64_t i = 0; i < nEntries; ++i ) {
            inTree->GetEntry(i);
            
            Cluster_map const *clusters
                = findClusterCenters( *strips, *cross_pos );

            for( int j = 0; j < event.nsim; ++j ) {
                if( event.simtype[j] == PHOTON
                    && event.simorig[j] == PION
                ) {
                    Photon photon(event, j);

                    // For this photon find the closest cluster overall
                    Cluster_map::const_iterator closestClusterIt
                        = clusters->cbegin();

                    // Start by assuming there is no cluster on each particular layer
                    // Find the closest by angDist, by phi, and by theta
                    Cluster_map::const_iterator closestClusterItLayer[LAYERS];
                    Cluster_map::const_iterator phiItLayer[LAYERS];
                    Cluster_map::const_iterator thetaItLayer[LAYERS];
                    for( int i = 0; i < LAYERS; ++i ) {
                        closestClusterItLayer[i] = clusters->cend();
                        phiItLayer[i] = clusters->cend();
                        thetaItLayer[i] = clusters->cend();
                    }

                    for( auto it = clusters->cbegin();
                         it != clusters->cend();
                         ++it
                    ) {
                        // Overall
                        if( angularDistance(cluster(it), photon) < angularDistance(cluster(closestClusterIt), photon) )
                            closestClusterIt = it;

                        // On the current layer
                        if( closestClusterItLayer[cluster(it).layer] == clusters->cend() )
                            closestClusterItLayer[cluster(it).layer] = it;
                        else if( angularDistance(cluster(it), photon) < angularDistance(cluster(closestClusterItLayer[cluster(it).layer]), photon ) )
                            closestClusterItLayer[cluster(it).layer] = it;

                        if( phiItLayer[cluster(it).layer] == clusters->cend() )
                            phiItLayer[cluster(it).layer] = it;
                        else if( std::fabs(cluster(it).cphi - photon.simphi) < std::fabs(cluster(phiItLayer[cluster(it).layer]).cphi - photon.simphi) )
                            phiItLayer[cluster(it).layer] = it;

                        if( thetaItLayer[cluster(it).layer] == clusters->cend() )
                            thetaItLayer[cluster(it).layer] = it;
                        else if( std::fabs(cluster(it).ctheta - photon.simtheta) < std::fabs(cluster(thetaItLayer[cluster(it).layer]).ctheta - photon.simtheta) )
                            thetaItLayer[cluster(it).layer] = it;
                    }

                    // If found no clusters whatsoever, give up on this photon
                    if( closestClusterIt == clusters->cend() ) {
                        continue;
                    }

                    // Fill the overall and layered histograms
                    double angDistValue = angularDistance( cluster(closestClusterIt), photon );
                    angDist->Fill( angDistValue );
                    
                    for( int i = 0; i < LAYERS; ++i ) {
                        if( closestClusterItLayer[i] == clusters->cend() )
                            continue;

                        double angDistValueLayer = angularDistance( cluster(closestClusterItLayer[i]), photon );
                        angDistLayer[i]->Fill( angDistValueLayer );
                    }
                    for( int i = 0; i < LAYERS; ++i ) {
                        if( phiItLayer[i] == clusters->cend() )
                            continue;

                        double phiDiffValueLayer = std::fabs( cluster(phiItLayer[i]).cphi - photon.simphi );
                        phiDistLayer[i]->Fill( phiDiffValueLayer );
                    }
                    for( int i = 0; i < LAYERS; ++i ) {
                        if( thetaItLayer[i] == clusters->cend() )
                            continue;

                        double thetaDiffValueLayer = std::fabs( cluster(thetaItLayer[i]).ctheta - photon.simtheta );
                        thetaDistLayer[i]->Fill( thetaDiffValueLayer );
                    }
                }
            }
        }
        
        std::cout << "\nAngular distance count: "
            << angDist->GetEntries() << std::endl;

        // Write histograms
        TFile *histoFile = new TFile("histos.root", "recreate");
        angDist->Write();
        for( int i = 0; i < LAYERS; ++i )
            angDistLayer[i]->Write();
        for( int i = 0; i < LAYERS; ++i )
            phiDistLayer[i]->Write();
        for( int i = 0; i < LAYERS; ++i )
            thetaDistLayer[i]->Write();
    }


    Cluster_map const *findClusterCenters(
        gera_nm::strip_data const &strips,
        gera_nm::cross_data const &cross_pos
    ) {
        Cluster_map *clusters = new Cluster_map();
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

            // Ensure that both clusters are on the same layer
            int layer1 = strips.layer[strip1];
            int layer2 = strips.layer[strip2];
            if (layer1 != layer2)
                continue;

            // Ensure that the current cluster pair is listed
            Cluster_map::iterator cl_it = clusters->find( cluster_id_paired );
            std::map< Cluster_id_t, Running >::iterator run_it = running.find( cluster_id_paired );
            if( cl_it == clusters->cend() ) {
                // If not found, add the pair to the list
                cl_it = ( clusters->insert( std::pair< Cluster_id_t, Cluster >(cluster_id_paired, Cluster()) ) ).first;
                cluster(cl_it).layer = layer1;
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
            cluster(cl_it).cphi
                = std::atan2(
                    run_it->second.y,
                    run_it->second.x
                );
            cluster(cl_it).ctheta
                = std::atan2(
                    std::sqrt(
                        std::pow(run_it->second.x, 2)
                        + std::pow(run_it->second.y, 2)
                    ),
                    run_it->second.z
                );
        }

        return clusters;
    }
    

    double angularDistance( Cluster const &cluster, Photon const &photon ) {
        double cphi = cluster.cphi;
        double ctheta = cluster.ctheta;
        double phphi = photon.simphi;
        double phtheta = photon.simtheta;

        double angle
            = std::acos(
                std::cos(ctheta)
                    * std::cos(phtheta)
                + std::sin(ctheta)
                    * std::sin(phtheta)
                    * std::cos(cphi - phphi)
            );
        return angle;
    }
}


template < typename T >
size_t findIndex(
    std::vector< T > const &vec,
    T const &elem
) {
    typename std::vector< T >::const_iterator it
        = std::find( vec.cbegin(), vec.cend(), elem );

    if( it == vec.cend() ) // not found
        return -1;
    return it - vec.cbegin();
}


void initLayeredHistos(
    TH1F **histArray,
    size_t size,
    std::string const &nameStart,
    std::string const &titleStart,
    size_t bins,
    double minX,
    double maxX
) {
    for( int i = 0; i < size; ++i ) {
        std::string name(nameStart), title(titleStart);
        name.push_back( '0' + i );
        title.push_back( '0' + i );
        histArray[i] = new TH1F( name.c_str(), title.c_str(), bins, minX, maxX );
    }
}


void selectPhotons(
    const std::string &inFileName = "/store25/semenov/strips_run039799.root",
    const std::string &outFileName = "/store25/bazhenov/piph.root"
) {
    cluster_div::selectPhotons( inFileName, outFileName );
}
