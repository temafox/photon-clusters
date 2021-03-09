#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>

#include <algorithm>
#include <iostream>

#include "sim_data.h"
#include "functions.h"
#include "cluster_div_config.h"

namespace cluster_div {

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
    auto *p_distanceFromClusterAxis = new MultipleHistos(2, "distanceFromClusterAxis", "Distance from cluster axis on layer: photons 1+", 25, 0, maxAngularDistance);
    auto *p_pionEnergy = new MultipleHistos(2, "pionEnergy", "Energy of the original pion", 20, PI_ENERGY_MIN, PI_ENERGY_MAX);
    auto &hitNumber = *p_hitNumber, &distanceFromClusterAxis = *p_distanceFromClusterAxis, &pionEnergy = *p_pionEnergy;

    // Process events
    //CUT_NENTRIES(1000);
    std::cout << "nEntries = " << nEntries << std::endl;
    for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
        inTree->GetEntry(entryIndex);

        ClusterMap *clusters = findClusterCenters(*strips, *cross_pos);

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

                inTree->GetEntry(cluster.eventIndex);
                pionEnergy[0].Fill(energy(PI_MC2, event.simmom[cluster.pionParticleIndex]));
            } else {
                distanceFromClusterAxis[1].Fill(std::max(angularDistance(cluster, cluster.photons[0]), angularDistance(cluster, cluster.photons[1])));
                hitNumber[1].Fill(cluster.hitNumber);

                inTree->GetEntry(cluster.eventIndex);
                pionEnergy[1].Fill(energy(PI_MC2, event.simmom[cluster.pionParticleIndex]));
            }
        }
    }

    std::cout << "\nDistance from cluster axis 1ph: " << distanceFromClusterAxis[0].GetEntries() << std::endl;
    std::cout << "Distance from cluster axis 2ph: " << distanceFromClusterAxis[1].GetEntries() << std::endl;

    std::cout << "\nHit number 1ph: " << hitNumber[0].GetEntries() << std::endl;
    std::cout << "Hit number 2ph: " << hitNumber[1].GetEntries() << std::endl;

    auto *c1 = new TCanvas();
    distanceFromClusterAxis[1].SetLineColor(kRed);
    distanceFromClusterAxis[0].Draw();
    distanceFromClusterAxis[1].Draw("same");
    c1->SetLogy();
    c1->Print("histos/distanceFromClusterAxis.png");

    auto *c2 = new TCanvas();
    hitNumber[1].SetLineColor(kRed);
    hitNumber[0].Draw();
    hitNumber[1].Draw("same");
    c2->SetLogy();
    c2->Print("histos/hitNumber.png");

    // Write histograms
    auto histoFile = new TFile(OUT_FILE_NAME, "recreate");
    for (int i = 0; i < 2; ++i) {
        hitNumber[i].Write();
        distanceFromClusterAxis[i].Write();
        pionEnergy[i].Write();
    }
}

}

void fillHistosForClassificationParameters(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    cluster_div::fillHistosForClassificationParameters(inFileName, outFileName);
}
