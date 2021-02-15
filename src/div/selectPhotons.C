#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

#include "sim_data.h"
#include "functions.h"

#define IN_FILE_NAME "/store25/semenov/strips_run039799.root"
#define OUT_FILE_NAME "/store25/bazhenov/piph.root"

namespace cluster_div {

void selectPhotons(
    const std::string &inFileName = IN_FILE_NAME,
    const std::string &outFileName = OUT_FILE_NAME
);

const ClusterMap *findClusterCenters(const gera_nm::strip_data &strips, const gera_nm::cross_data &cross_pos);

}
    
/// Implementations

namespace cluster_div {

void selectPhotons(const std::string &inFileName, const std::string &outFileName) {
    // Initialize inputs
    auto inFile = new TFile(inFileName.c_str());
    auto inTree = (TTree *)inFile->Get("tr_lxe;80");

    Long64_t nEntries = inTree->GetEntries();

    gera_nm::tree_data event{};
    auto strips = new gera_nm::strip_data;
    auto cross_pos = new gera_nm::cross_data;

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
    TH1F *angDist = new TH1F("angDist", "Angular distance to closest cluster", 70, 0, TMath::Pi());

    LayeredHistos angleLayeredHistos(LAYERS, "angDist", "Angular distance at layer ", 70, 0, TMath::Pi());
    LayeredHistos phiLayeredHistos(LAYERS, "phiDist", "Difference for phi at layer ", 140, 0, 2. * TMath::Pi());
    LayeredHistos thetaLayeredHistos(LAYERS, "thetaDist", "Difference for theta at layer ", 70, 0, TMath::Pi());

    // Angular distributions of pions
    TH1F *piTheta = new TH1F("piTheta", "pi_simtheta", 70, 0, TMath::Pi());
    TH1F *piPhi = new TH1F("piPhi", "pi_simphi", 140, 0, 2. * TMath::Pi());

    // Process events
    nEntries = 100;
    std::cout << "nEntries = " << nEntries << std::endl;
    for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
        inTree->GetEntry(entryIndex);

        const ClusterMap *clusters = findClusterCenters(*strips, *cross_pos);

        for (int particleIndex = 0; particleIndex < event.nsim; ++particleIndex) {
            // pi0 -> _2gamma_
            if ((event.simtype[particleIndex] == PHOTON) && (event.simorig[particleIndex] == PION)) {
                Photon photon(event, particleIndex);
                auto closestClusterIt = clusters->cbegin();

                // Start by assuming there is no cluster on each particular layer
                // Find the closest by angDist, by phi, and by theta
                ClusterMap::const_iterator layeredClosestClusterIt[LAYERS];
                ClusterMap::const_iterator layeredClosestPhiIt[LAYERS];
                ClusterMap::const_iterator layeredClosestThetaIt[LAYERS];
                for (int k = 0; k < LAYERS; ++k ) {
                    layeredClosestClusterIt[k] = clusters->cend();
                    layeredClosestPhiIt[k] = clusters->cend();
                    layeredClosestThetaIt[k] = clusters->cend();
                }

                for (auto it = clusters->cbegin(); it != clusters->cend(); ++it) {
                    auto &cluster = it->second;
                    auto &closestCluster = closestClusterIt->second;
                    auto &layeredClosestCluster = layeredClosestClusterIt[cluster.layer]->second;
                    auto &layeredClosestPhi = layeredClosestPhiIt[cluster.layer]->second;
                    auto &layeredClosestTheta = layeredClosestThetaIt[cluster.layer]->second;
                    
                    // Overall
                    if (angularDistance(cluster, photon) < angularDistance(closestCluster, photon))
                        closestClusterIt = it;

                    // On the current layer
                    if (layeredClosestClusterIt[cluster.layer] == clusters->cend())
                        layeredClosestClusterIt[cluster.layer] = it;
                    else if (angularDistance(cluster, photon) < angularDistance(layeredClosestCluster, photon))
                        layeredClosestClusterIt[cluster.layer] = it;

                    if (layeredClosestPhiIt[cluster.layer] == clusters->cend())
                        layeredClosestPhiIt[cluster.layer] = it;
                    else if (std::fabs(cluster.cphi - photon.simphi) < std::fabs(layeredClosestPhi.cphi - photon.simphi))
                        layeredClosestPhiIt[cluster.layer] = it;

                    if (layeredClosestThetaIt[cluster.layer] == clusters->cend())
                        layeredClosestThetaIt[cluster.layer] = it;
                    else if (std::fabs(cluster.ctheta - photon.simtheta) < std::fabs(layeredClosestTheta.ctheta - photon.simtheta))
                        layeredClosestThetaIt[cluster.layer] = it;
                }

                // If we found no clusters, give up on this photon
                if (closestClusterIt == clusters->cend())
                    continue;

                // Fill the overall and layered histograms
                double angDistValue = angularDistance(closestClusterIt->second, photon);
                angDist->Fill(angDistValue);

                for (int k = 0; k < LAYERS; ++k) {
                    if (layeredClosestClusterIt[k] == clusters->cend())
                        continue;

                    double angDistValueLayer = angularDistance(layeredClosestClusterIt[k]->second, photon);
                    angleLayeredHistos[k].Fill(angDistValueLayer);
                }
                for (int k = 0; k < LAYERS; ++k) {
                    if (layeredClosestPhiIt[k] == clusters->cend())
                        continue;

                    double phiDiffValueLayer = std::fabs(layeredClosestPhiIt[k]->second.cphi - photon.simphi);
                    phiLayeredHistos[k].Fill(phiDiffValueLayer);
                }
                for (int k = 0; k < LAYERS; ++k) {
                    if (layeredClosestThetaIt[k] == clusters->cend())
                        continue;

                    double thetaDiffValueLayer = std::fabs(layeredClosestThetaIt[k]->second.ctheta - photon.simtheta);
                    thetaLayeredHistos[k].Fill(thetaDiffValueLayer);
                }
            }

            // e+ e- -> _pi0_ gamma
            if (event.simtype[particleIndex] == PION) {
                piTheta->Fill(event.simtheta[particleIndex]);
                piPhi->Fill(event.simphi[particleIndex]);
            }


        }
    }

    std::cout << "\nAngular distance count: " << angDist->GetEntries() << std::endl;

    // Write histograms
    auto histoFile = new TFile("histos.root", "recreate");
    angDist->Write();
    for (int i = 0; i < LAYERS; ++i)
        angleLayeredHistos[i].Write();
    for (int i = 0; i < LAYERS; ++i)
        phiLayeredHistos[i].Write();
    for (int i = 0; i < LAYERS; ++i)
        thetaLayeredHistos[i].Write();
    piTheta->Write();
    piPhi->Write();
}


const ClusterMap *findClusterCenters(const gera_nm::strip_data &strips, const gera_nm::cross_data &cross_pos) {
    auto *clusters = new ClusterMap();
    size_t numberOfCrosses = cross_pos.id1.size();

    class Running {
    public:
        double x, y, z;
        size_t size;

        Running(): x(0.), y(0.), z(0.), size(0) {}
    };
    std::map<Cluster_id_t, Running> running;

    for (size_t cross_i = 0; cross_i < numberOfCrosses; ++cross_i) {
        // Find the indices of the crossing strips
        size_t strip1 = findIndex(strips.packedID, cross_pos.id1[cross_i]);
        size_t strip2 = findIndex(strips.packedID, cross_pos.id2[cross_i]);

        // Ensure that either cluster_id != -1 (means no cluster)
        int cluster_id1 = strips.cluster_id[strip1];
        int cluster_id2 = strips.cluster_id[strip2];
        if ((cluster_id1 == -1) || (cluster_id2 == -1))
            continue;
        Cluster_id_t cluster_id_paired = Cluster_id_t(cluster_id1, cluster_id2);

        // Ensure that both clusters are on the same layer
        int layer1 = strips.layer[strip1];
        int layer2 = strips.layer[strip2];
        if (layer1 != layer2)
            continue;

        // Ensure that the current cluster pair is listed
        auto cl_it = clusters->find(cluster_id_paired);
        auto run_it = running.find(cluster_id_paired);
        if (cl_it == clusters->cend()) {
            // If not found, add the pair to the list
            cl_it = (clusters->insert(std::pair<Cluster_id_t, Cluster>(cluster_id_paired, Cluster()))).first;
            cl_it->second.layer = layer1;
            run_it = (running.insert(std::pair<Cluster_id_t, Running>(cluster_id_paired, Running()))).first;
        }

        run_it->second.x += cross_pos.x[cross_i];
        run_it->second.y += cross_pos.y[cross_i];
        run_it->second.z += cross_pos.z[cross_i];
        run_it->second.size += 1;
    }

    auto cl_it = clusters->begin();
    auto run_it = running.begin();
    for (; cl_it != clusters->end() && run_it != running.cend(); ++cl_it, ++run_it) {
        run_it->second.x /= run_it->second.size;
        run_it->second.y /= run_it->second.size;
        run_it->second.z /= run_it->second.size;

        // tg(phi) = y/x
        // tg(theta) = r/z, r^2 = x^2 + y^2
        cl_it->second.cphi = std::atan2(run_it->second.y,run_it->second.x);
        cl_it->second.ctheta = std::atan2(std::sqrt(std::pow(run_it->second.x, 2) + std::pow(run_it->second.y, 2)), run_it->second.z);
    }

    return clusters;
}

}

void selectPhotons(const std::string &inFileName, const std::string &outFileName) {
    cluster_div::selectPhotons( inFileName, outFileName );
}
