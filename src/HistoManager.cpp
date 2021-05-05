#include <stdexcept>
#include "HistoManager.h"

void HistoManager::addHisto(const std::string &name, std::size_t bins, double xmin, double xmax) {
    auto histo = new TH1F(name.c_str(), name.c_str(), bins, xmin, xmax);
    singleHistos_.insert({name, histo});
}


std::string buildSerialName(const std::string &namePrefix, std::size_t index, std::size_t arraySize) {
    int digitCount = 1;
    while (arraySize > 0) {
        arraySize /= 10;
        digitCount++;
    }

    std::string result = namePrefix;
    result.resize(result.size() + digitCount);
    for (auto i = 0; i < digitCount; ++i) {
        result.append();
    }
}

void HistoManager::addArray(std::size_t size, const std::string &namePrefix, std::size_t bins, double xmin,
                            double xmax) {
    multipleHistos_.insert({namePrefix, std::vector<TH1F *>(size)});
    auto &histoArray = multipleHistos_[namePrefix];

    for (std::size_t i = 0; i < size; ++i) {

        histoArray[i] = new TH1F();
    }
}


TH1F *HistoManager::getHisto(const std::string &name) {
    return singleHistos_[name];
}


TH1F *HistoManager::getHisto(const std::string &namePrefix, std::size_t index) {

}