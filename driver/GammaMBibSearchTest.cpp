#include <fstream>
#include <NuclidLibrary.h>
#include <BibSearchCtx.h>
#include <assert.h>
#include "../finding/GammaMBibSearch.h"
#include <cmath>
#include "../baseline/BinBaseline.h"

class GammaMBibSearchTest{
    public:

    void initialize_testing_environment();
    double get_median_diff_from_expected(const FitOut& fit, const Nuclid& nucl);
    bool test_gamma_m_bib_search();

    private:
    BibSearchCtx ctx;
};

void GammaMBibSearchTest::initialize_testing_environment(){

    std::vector<double> counts;
    std::ifstream nuclideFile("path/to/spec.txt");
    if (!nuclideFile.is_open()) {
        throw std::runtime_error("Could not open nuclide library file.");
    }
    double energy;
    while(nuclideFile >> energy){
        counts.push_back(energy);
    }

    NuclidLibrary nuclideLibrary;
    nuclideLibrary.loadLibrary("path/to/nuclide_library.txt");

    FWHMC fwhm_model = FWHMC{1.0, 0.01, 1e-6};
    EC ec_model = EC{0.0, 0.5, 1e-4};
    std::vector<Nuclid> expected_nuclides = nuclideLibrary.getNuclideList();

    std::vector<double> baseline = BinBaseline::estimate(ctx.counts);

    ctx = BibSearchCtx{expected_nuclides, fwhm_model, ec_model, counts, baseline, 1e-3, 8};

}

bool GammaMBibSearchTest::test_gamma_m_bib_search(){

    assert(!ctx.expected_nuclides.empty());
    assert(ctx.Q0_conv_tol > 0);
    assert(ctx.iter > 0);
    assert(!ctx.counts.empty());


    std::vector<std::pair<Nuclid, FitOut>> results = GammaMBibSearch::search(ctx);
    assert(!results.empty());

    std::ofstream outputFile("path/to/output.txt");
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file.");
    }
    for (const auto& fit : results) {
        outputFile << "Nuclide: Z=" << fit.first.atomic_number << " A=" << fit.first.mass_number << "\n";
        outputFile << "B0:" << fit.second.b0 << "  B1:" << fit.second.b1 << " Q0:" << fit.second.Q0 << " Counts(Center):" << fit.second.counts_center << "\n";
        outputFile << "(Chi^2)/dof: " << fit.second.reduced_chi2 << "\n";
        outputFile << "-----------------------------------\n";

    }


    return true;
}
