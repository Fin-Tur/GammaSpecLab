//
// Created by f.willems on 18.11.2025.
//

#include "SecondDSearchTest.h"

#include <stdexcept>
#include <vector>
#include <fstream>
#include <iostream>

#include "../models/cal/EC.h"
#include "../models/cal/FWHMC.h"
#include "../models/elements/Peak.h"
#include "../models/ctx/SecDervCtx.h"
#include "../finding/SecondDerivateSearch.h"

void SecondDSearchTest::test_seconddsearch() {
    std::vector<double> counts;
    std::ifstream specFile("C:/Users/f.willems/CLionProjects/Toolbox/data/Co-2106.TKA");
    if (!specFile.is_open()) {
        throw std::runtime_error("Could not open spectrum file.");
    }

    double count;
    while (specFile >> count) {
        counts.push_back(count);
    }
    specFile.close();

    auto fwhm_model = FWHMC{3.930e-001, 3.523e-002, 0.0};
    auto ec_model = EC{-2.694e-001, 1.982e-001, 0.0};

    const SecDervCtx ctx{counts, fwhm_model, ec_model, 3, 100};

    std::vector<Peak> found_peaks = SecondDerivateSearch::find_peaks(ctx);

    for (const auto& peak : found_peaks) {
        std::cout << "Found peak at " << peak.peak_center << " with peak height: " << peak.peak_height << " and z: " << peak.z << std::endl;
    }
    return;

}
