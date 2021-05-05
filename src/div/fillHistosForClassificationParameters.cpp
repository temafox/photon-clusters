#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>

#include <algorithm>
#include <iostream>
#include <list>

#include "sim_data.h"
#include "functions.h"
#include "cluster_div_config.h"

namespace cluster_div {

void drawAndSaveTwoHistos(MultipleHistos &histoArray, const char *fileName) {
    auto *c = new TCanvas();

    histoArray[0].Scale(1./histoArray[0].GetEntries());
    histoArray[1].Scale(1./histoArray[1].GetEntries());
    histoArray[1].SetLineColor(kRed);
    histoArray[0].Draw();
    histoArray[1].Draw("same");

    c->SetLogy();

    c->Print(fileName);
}

std::list<Cluster_id_t> findClusterChains(const ClusterMap &clusters, double maxAngularDistance) {
    std::list<Cluster_id_t> res;

    for (const auto it : clusters) {
        if (!res.empty()) {
            if ((it.second.layer > clusters[res.back()].layer) &&
                    (angularDistance(it.second, clusters[res.back()]) < maxAngularDistance))
                res.push_back(it.first);
        } else {
            ;
        }
    }

    return res;
}

void fillHistosForClassificationParameters(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
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

    // Parameters
    const double maxAngularDistance = TMath::Pi();
    auto *p_hitNumber = new MultipleHistos(2, "hitNumber", "Number of hits on layer: photons 1+", 25, 0, 300);
    auto &hitNumber = *p_hitNumber;

    auto *p_distanceFromClusterAxis = new MultipleHistos(2, "distanceFromClusterAxis", "Distance from cluster axis on layer: photons 1+", 25, 0, maxAngularDistance);
    auto &distanceFromClusterAxis = *p_distanceFromClusterAxis;

    auto *p_clusterArea = new MultipleHistos(2, "clusterArea", "Area of cluster on layer 6: photons 1+", 20, 0, 100);
    auto &clusterArea = *p_clusterArea;

    auto *p_clusterCenterDisplacement = new MultipleHistos(2, "clusterCenterDisplacement", "Displacement of cluster center relative to previous layer: photons 1+", 20, 0, 100);
    auto &clusterCenterDisplacement = *p_clusterCenterDisplacement;

    auto *p_pionEnergy = new MultipleHistos(2, "pionEnergy", "Energy of the original pion: photons 1+", 20, PI_ENERGY_MIN, PI_ENERGY_MAX);
    auto &pionEnergy = *p_pionEnergy;

    // Process events
    //CUT_NENTRIES(1000);
    std::cout << "nEntries = " << nEntries << std::endl;
    for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
        inTree->GetEntry(entryIndex);

        ClusterMap *clusters = findClusterCenters(*strips, *cross_pos);
        std::list<Cluster_id_t> clusterChains = findClusterChains(*clusters, maxAngularDistance);

        for (int particleIndex = 0; particleIndex < event.nsim; ++particleIndex) {
            // pi0 -> _2gamma_
            if ((event.simtype[particleIndex] == PHOTON) && (event.simorig[particleIndex] == PION)) {
                auto closestClusterIt = clusters->begin();
                if (closestClusterIt == clusters->cend())
                    continue;

                Photon photon(event, particleIndex);

                for (auto it = clusters->begin(); it != clusters->end(); ++it) {
                    auto &cluster = it->second;
                    auto &closestCluster = closestClusterIt->second;

                    // Overall
                    if (angularDistance(cluster, photon) < angularDistance(closestCluster, photon))
                        closestClusterIt = it;
                }

                if (angularDistance(closestClusterIt->second, photon) < maxAngularDistance) {
                    closestClusterIt->second.photons.push_back(photon);
                    closestClusterIt->second.eventIndex = entryIndex;
                    closestClusterIt->second.pionParticleIndex = findIndex(event.simtype, event.nsim, PION);
                }
            }
        }

        for (const auto &clusterKvp : *clusters) {
            auto &cluster = clusterKvp.second;
            
            if (cluster.photons.empty())
                continue;
            else if (cluster.photons.size() == 1) {
                distanceFromClusterAxis[0].Fill(angularDistance(cluster, cluster.photons[0]));
                hitNumber[0].Fill(cluster.hitNumber);

                if (cluster.layer > 0) {
                    clusterCenterDisplacement[0].Fill(0);//angularDistance(cluster, prevCluster));
                }

                inTree->GetEntry(cluster.eventIndex);
                pionEnergy[0].Fill(energy(PI_MC2, event.simmom[cluster.pionParticleIndex]));

                clusterArea[0].Fill(0);
            } else {
                distanceFromClusterAxis[1].Fill(std::max(angularDistance(cluster, cluster.photons[0]), angularDistance(cluster, cluster.photons[1])));
                hitNumber[1].Fill(cluster.hitNumber);

                if (cluster.layer > 0) {
                    clusterCenterDisplacement[1].Fill(0);//angularDistance(cluster, prevCluster));
                }

                inTree->GetEntry(cluster.eventIndex);
                pionEnergy[1].Fill(energy(PI_MC2, event.simmom[cluster.pionParticleIndex]));

                clusterArea[1].Fill(0);
            }
        }
    }

    std::cout << "\nDistance from cluster axis 1ph: " << distanceFromClusterAxis[0].GetEntries() << std::endl;
    std::cout << "Distance from cluster axis 2ph: " << distanceFromClusterAxis[1].GetEntries() << std::endl;

    std::cout << "\nHit number 1ph: " << hitNumber[0].GetEntries() << std::endl;
    std::cout << "Hit number 2ph: " << hitNumber[1].GetEntries() << std::endl;

    drawAndSaveTwoHistos(distanceFromClusterAxis, "histos/distanceFromClusterAxis.png");
    drawAndSaveTwoHistos(hitNumber, "histos/hitNumber.png");
    drawAndSaveTwoHistos(clusterCenterDisplacement, "histos/clusterCenterDisplacement.png");
    drawAndSaveTwoHistos(clusterArea, "histos/clusterArea.png");

    // Write histograms
    auto histoFile = new TFile(OUT_FILE_NAME, "recreate");
    for (int i = 0; i < 2; ++i) {
        hitNumber[i].Write();
        distanceFromClusterAxis[i].Write();
        clusterCenterDisplacement[i].Write();
        clusterArea[i].Write();
        pionEnergy[i].Write();
    }
}

}

void fillHistosForClassificationParameters(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    cluster_div::fillHistosForClassificationParameters(inFileName, outFileName);
}
