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

#define IN_FILE_NAME "/store25/semenov/strips_run039799.root"
#define OUT_FILE_NAME "class_params.root"

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
    MultipleHistos hitNumber(2, "hitNumber", "Number of hits on layer: photons 1+", 70, 0, 100);
    MultipleHistos distanceFromClusterAxis(2, "distanceFromClusterAxis", "Distance from cluster axis on layer: photons 1+", 70, 0, TMath::Pi());
    const double maxAngularDistance = 0.3;

    // Process events
    CUT_NENTRIES;
    std::cout << "nEntries = " << nEntries << std::endl;
    for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
        inTree->GetEntry(entryIndex);

        ClusterMap *clusters = findClusterCenters(*strips, *cross_pos);

        for (int particleIndex = 0; particleIndex < event.nsim; ++particleIndex) {
            // pi0 -> _2gamma_
            if ((event.simtype[particleIndex] == PHOTON) && (event.simorig[particleIndex] == PION)) {
                Photon photon(event, particleIndex);
                auto closestClusterIt = clusters->begin();

                for (auto it = clusters->begin(); it != clusters->end(); ++it) {
                    auto &cluster = it->second;
                    auto &closestCluster = closestClusterIt->second;

                    // Overall
                    if (angularDistance(cluster, photon) < angularDistance(closestCluster, photon))
                        closestClusterIt = it;
                }

                if (angularDistance(closestClusterIt->second, photon) < maxAngularDistance)
                    closestClusterIt->second.photons.push_back(photon);
            }
        }

        for (auto clusterIt = clusters->cbegin(); clusterIt != clusters->cend(); ++clusterIt) {
            if (clusterIt->second.photons.size() == 1)
                distanceFromClusterAxis[0].Fill(angularDistance(clusterIt->second, clusterIt->second.photons[0]));
            else
                distanceFromClusterAxis[1].Fill(std::max(angularDistance(clusterIt->second, clusterIt->second.photons[0]), angularDistance(clusterIt->second, clusterIt->second.photons[1])));
        }
    }

    distanceFromClusterAxis[0].Draw();
    distanceFromClusterAxis[1].Draw("same");
}

}

void fillHistosForClassificationParameters(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    cluster_div::fillHistosForClassificationParameters(inFileName, outFileName);
}
