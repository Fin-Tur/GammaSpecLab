#include <fstream>
#include <NuclidLibrary.h>
#include <BibSearchCtx.h>
#include <assert.h>
#include "../finding/GammaMBibSearch.h"
#include <cmath>

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

    ctx = BibSearchCtx{expected_nuclides, fwhm_model, ec_model, counts, 1e-3, 8};

}

double GammaMBibSearchTest::get_median_diff_from_expected(const FitOut& fit, const Nuclid& nucl){
        int ck = ctx.ec_model.channel_at(nucl.gamma_energy);
        double Ec = ctx.ec_model.energy_at(ck);
        double FWHM_keV = ctx.fwhm_model.fwhm_at(Ec);
        double z_keV = FWHM_keV / (std::sqrt(2.0*std::log(2)));

    for(int i = ctx.bg_ctx.lR; i <= ctx.bg_ctx.rR; ++i){

    }
    return 0.0;
}

bool GammaMBibSearchTest::test_gamma_m_bib_search(){

    assert(!ctx.expected_nuclides.empty());
    assert(ctx.Q0_conv_tol > 0);
    assert(ctx.iter > 0);
    assert(!ctx.counts.empty());

    std::vector<FitOut> results = GammaMBibSearch::search(ctx);
    assert(!results.empty());

    std::ofstream outputFile("path/to/output.txt");
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open output file.");
    }
    for (const auto& fit : results) {
        outputFile << "B0:" << fit.b0 << "  B1:" << fit.b1 << " Q0:" << fit.Q0 << " Height:" << fit.counts_center << "\n";   
        outputFile << "Median difference from expected: " << get_median_diff_from_expected(fit) << "\n"; //FIXME need nuclide here?
        outputFile << "-----------------------------------\n";

    }


    return true;
}
