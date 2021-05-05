#include <string>
#include "sim_data.h"
#include "cluster_div_config.h"
#include "ClusterParser.h"
#include "HistoManager.h"

void classParams(const std::string &inFileName = IN_FILE_NAME, const std::string &outFileName = OUT_FILE_NAME) {
    auto clusterParser = ClusterParser{inFileName};

    auto histoManager = HistoManager{};
    histoManager.addArray(2, "hitNumber", 25, 0., 300.);
    histoManager.addArray(2, "clusterArea", 20, 0., 100.);
    histoManager.addArray(2, "clusterCenterDisplacement", 20, 0., 100.);
    histoManager.addArray(2, "deviationRatio", 20, 0., 10.);

    for (auto event : clusterParser) {
        for (auto particle : event) {
            if ((particle.getType() == PHOTON) && (particle.getOrigin() == PION)) {
                ;
            }
        }
    }
}