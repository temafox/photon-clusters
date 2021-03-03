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
#include "cluster_div_config.h"

namespace cluster_div {

void selectPhotons(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    // Initialize inputs
    auto inFile = new TFile(inFileName.c_str());
    auto inTree = (TTree *)inFile->Get("tr_lxe");

    Long64_t nEntries = inTree->GetEntries();

    gera_nm::tree_data event{};
    auto strips = new gera_nm::strip_data;
    auto cross_pos = new gera_nm::cross_data;

    inTree->SetBranchAddress("nsim", &(event.nsim));
    inTree->SetBranchAddress("simtype", event.simtype);
    inTree->SetBranchAddress("simorig", event.simorig);
    inTree->SetBranchAddress("simmom", event.simmom);
    inTree->SetBranchAddress("simphi", event.simphi);
    inTree->SetBranchAddress("simtheta", event.simtheta);
    inTree->SetBranchAddress("simvtx", event.simvtx);
    inTree->SetBranchAddress("simvty", event.simvty);
    inTree->SetBranchAddress("simvtz", event.simvtz);
    inTree->SetBranchAddress("strips", &strips);
    inTree->SetBranchAddress("cross_pos", &cross_pos);

    // Create histograms of angular distances
    // Layers are numbered 0..6
    TH1F *angDist = new TH1F("angDist", "Angular distance to closest cluster", 70, 0, TMath::Pi());

    MultipleHistos angleLayeredHistos(LAYERS, "angDist", "Angular distance at layer ", 70, 0, TMath::Pi());
    MultipleHistos phiLayeredHistos(LAYERS, "phiDist", "Difference for phi at layer ", 140, 0, TMath::Pi());
    MultipleHistos thetaLayeredHistos(LAYERS, "thetaDist", "Difference for theta at layer ", 70, 0, TMath::Pi());

    // Angular and energy distributions of pions
    TH1F *piTheta = new TH1F("piTheta", "pi_simtheta", 70, 0, TMath::Pi());
    TH1F *piPhi = new TH1F("piPhi", "pi_simphi", 140, 0, 2. * TMath::Pi());
    TH1F *piEnergy = new TH1F("piEnergy", "pi_energy", 70, PI_ENERGY_MIN, PI_ENERGY_MAX);

    // Process events
    //CUT_NENTRIES(100);
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
                    else if (phiDistance(cluster.cphi, photon.simphi) < phiDistance(layeredClosestPhi.cphi, photon.simphi))
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

                    double phiDiffValueLayer = phiDistance(layeredClosestPhiIt[k]->second.cphi, photon.simphi);
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
                piEnergy->Fill(energy(PI_MC2, event.simmom[particleIndex]));
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
    piEnergy->Write();
}

}

void selectPhotons(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    cluster_div::selectPhotons(inFileName, outFileName);
}
