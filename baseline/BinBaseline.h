//
// Created by f.willems on 13.10.2025.
//

#ifndef TOOLBOX_BINBASELINE_H
#define TOOLBOX_BINBASELINE_H
#include <vector>
#include"..//models/Bin.h"

namespace BinBaseline {
    std::vector<double> estimate(const std::vector<double> &counts);
};


#endif //TOOLBOX_BINBASELINE_H