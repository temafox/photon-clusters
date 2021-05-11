#ifndef CLUSTER_DIV_MULTIPLEHISTOS_H
#define CLUSTER_DIV_MULTIPLEHISTOS_H

#include <TH1F.h>

#include <string>

namespace cluster_div {

class MultipleHistos {
public:
    TH1F **histArray;
    size_t size;
    std::string nameStart;
    std::string titleStart;
    size_t bins;
    double minX;
    double maxX;

    MultipleHistos(size_t size, std::string nameStart, std::string titleStart, size_t bins, double minX, double maxX);
    ~MultipleHistos();

    TH1F &operator[](size_t layer) const;
private:
    void initMultipleHistos();
};

/// Implementation

MultipleHistos::MultipleHistos(size_t size, std::string nameStart, std::string titleStart, size_t bins, double minX, double maxX):
        size(size),
        nameStart(std::move(nameStart)),
        titleStart(std::move(titleStart)),
        bins(bins),
        minX(minX),
        maxX(maxX)
{
    histArray = new TH1F*[size];
    initMultipleHistos();
}

MultipleHistos::~MultipleHistos() {
    for (int i = 0; i < size; ++i)
        delete histArray[i];
    delete[] histArray;
}

TH1F &MultipleHistos::operator[](size_t layer) const {
    return *(histArray[layer]);
}

void MultipleHistos::initMultipleHistos() {
    histArray = new TH1F*[size];

    for (int i = 0; i < size; ++i) {
        std::string name(nameStart), title(titleStart);
        name.append(std::to_string(i));
        title.append(std::to_string(i));

        histArray[i] = new TH1F(name.c_str(), title.c_str(), bins, minX, maxX);
    }
}

}

#endif // CLUSTER_DIV_MULTIPLEHISTOS_H
