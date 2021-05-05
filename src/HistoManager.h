#ifndef CLUSTER_DIV_HISTOMANAGER_H
#define CLUSTER_DIV_HISTOMANAGER_H

#include <string>
#include <map>
#include <vector>
#include <TH1F.h>

// Memory leaks are intentional!!!
// This is due to ROOT's peculiarities
class HistoManager {
private:
    std::map<std::string, TH1F *> singleHistos_{};
    std::map<std::string, std::vector<TH1F *>> multipleHistos_{};

public:
    void addHisto(const std::string &name, std::size_t bins, double xmin, double xmax);
    void addArray(std::size_t size, const std::string &namePrefix, std::size_t bins, double xmin, double xmax);

    TH1F *getHisto(const std::string &name);
    TH1F *getHisto(const std::string &namePrefix, std::size_t index);
};


#endif //CLUSTER_DIV_HISTOMANAGER_H
