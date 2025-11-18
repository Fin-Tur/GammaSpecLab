#include "GammaMBibSearchTest.h"

#include <fstream>
#include <cassert>
#include <cmath>
#include <iostream>

#include "../baseline/AlsBaseline.h"
#include "../baseline/BinBaseline.h"
#include "../nuclides/NuclidLibrary.h"
#include "../finding/GammaMBibSearch.h"



GammaMBibSearchTest::GammaMBibSearchTest() : ctx{} {
    //this->ctx = BibSearchCtx{};
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

    NuclidLibrary nuclideLibrary;
    nuclideLibrary.loadLibrary("C:/Users/f.willems/CLionProjects/Toolbox/data/nuclide_library.txt");

    auto fwhm_model = FWHMC{3.930e-001, 3.523e-002, 0.0};
    auto ec_model = EC{-2.694e-001, 1.982e-001, 0.0};

    std::vector<Nuclid> expected_nuclides = nuclideLibrary.getNuclideList();

    std::vector<double> baseline = BinBaseline::estimate(counts, fwhm_model);
    //std::vector<double> baseline = ALS_BASELINE::estimateBackgroundUsingALS();

    this->ctx = BibSearchCtx{expected_nuclides, fwhm_model, ec_model, counts, baseline, 1e-3, 8};
}


bool GammaMBibSearchTest::test_gamma_m_bib_search(){

    assert(!ctx.expected_nuclides.empty());
    assert(ctx.Q0_conv_tol > 0);
    assert(ctx.iter > 0);
    assert(!ctx.counts.empty());


    std::vector<std::pair<Nuclid, FitOut>> results = GammaMBibSearch::search(ctx);
    assert(!results.empty());

    std::ofstream outputFile("test.txt");
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file.");
    }
    
    for (const auto& fit : results) {
        outputFile << "Nuclide: Z=" << fit.first.atomic_number << " A=" << fit.first.mass_number << "\n";
        outputFile << "B0:" << fit.second.b0 << "  B1:" << fit.second.b1 << " Q0:" << fit.second.Q0 << " Counts(Center):" << fit.second.counts_center << "\n";
        outputFile << "(Chi^2)/dof: " << fit.second.reduced_chi2 << "\n";
        outputFile << "-----------------------------------\n";

    }

    outputFile.close();
    return true;
}


