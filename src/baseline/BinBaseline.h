//
// Created by f.willems on 13.10.2025.
//

#pragma once

#include <vector>

#include "../models/cal/FWHMC.h"

namespace BinBaseline {
    std::vector<double> estimate(const std::vector<double> &counts, FWHMC fwhm_model);
};

