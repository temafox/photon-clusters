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
    const double maxAngularDistance = 0.5;
    MultipleHistos *p_hitNumber = new MultipleHistos(2, "hitNumber", "Number of hits on layer: photons 1+", 25, 0, 100);
    MultipleHistos *p_distanceFromClusterAxis = new MultipleHistos(2, "distanceFromClusterAxis", "Distance from cluster axis on layer: photons 1+", 25, 0, maxAngularDistance);
    MultipleHistos &hitNumber = *p_hitNumber, &distanceFromClusterAxis = *p_distanceFromClusterAxis;

    // Process events
    CUT_NENTRIES(1000);
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

                if (angularDistance(closestClusterIt->second, photon) < maxAngularDistance)
                    closestClusterIt->second.photons.push_back(photon);
            }
        }

        for (auto clusterIt = clusters->cbegin(); clusterIt != clusters->cend(); ++clusterIt) {
            if (clusterIt->second.photons.size() == 0)
                continue;
            else if (clusterIt->second.photons.size() == 1) {
                distanceFromClusterAxis[0].Fill(angularDistance(clusterIt->second, clusterIt->second.photons[0]));
                hitNumber[0].Fill(clusterIt->second.hitNumber);
            } else {
                distanceFromClusterAxis[1].Fill(std::max(angularDistance(clusterIt->second, clusterIt->second.photons[0]), angularDistance(clusterIt->second, clusterIt->second.photons[1])));
                hitNumber[1].Fill(clusterIt->second.hitNumber);
            }
        }
    }

    std::cout << "Distance from cluster axis 1ph: " << distanceFromClusterAxis[0].GetEntries() << std::endl;
    std::cout << "Distance from cluster axis 2ph: " << distanceFromClusterAxis[1].GetEntries() << std::endl;

    std::cout << "Hit number 1ph: " << hitNumber[0].GetEntries() << std::endl;
    std::cout << "Hit number 2ph: " << hitNumber[1].GetEntries() << std::endl;

    TCanvas *c1 = new TCanvas();
    distanceFromClusterAxis[1].SetLineColor(kRed);
    distanceFromClusterAxis[0].Draw();
    distanceFromClusterAxis[1].Draw("same");
    gPad->SetLogy();
    c1->Print("dfca.png");

    TCanvas *c2 = new TCanvas();
    hitNumber[1].SetLineColor(kRed);
    hitNumber[0].Draw();
    hitNumber[1].Draw("same");
    gPad->SetLogy();
    c2->Print("hn.png");
}

}

void fillHistosForClassificationParameters(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    cluster_div::fillHistosForClassificationParameters(inFileName, outFileName);
}
