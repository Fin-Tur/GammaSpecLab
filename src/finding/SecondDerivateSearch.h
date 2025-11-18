#pragma once
#include <vector>
#include "../models/ctx/SecDervCtx.h"
#include "../models/elements/Peak.h"

namespace SecondDerivateSearch {

    std::vector<Peak> find_peaks(const SecDervCtx& ctx);
    Peak construct_peak(const SecDervCtx& ctx, const std::vector<std::pair<int, double>>& ss_vals);
    
};