#include "../models/ctx/BibSearchCtx.h"
#include <fstream>
#include "../nuclides/NuclidLibrary.h"
#include "../baseline/BinBaseline.h"
#include "../finding/GammaMBibSearch.h"

class GammaMBibSearchDriver{

    public:
        void initialize_driver_environment(const std::string& spec_path, const std::string& nuclide_lib_path, const FWHMC& fwhm_model, const EC& ec_model, double Q0_conv_tol, int iter);
        std::vector<Nuclid> run(const double chi2_threshold);

    private:
        BibSearchCtx ctx;
 };

 void GammaMBibSearchDriver::initialize_driver_environment(const std::string& spec_path, const std::string& nuclide_lib_path, const FWHMC& fwhm_model, const EC& ec_model, const double Q0_conv_tol, const int iter){

    std::vector<double> counts;
    std::ifstream nuclideFile(spec_path);
    if (!nuclideFile.is_open()) {
        throw std::runtime_error("Could not open nuclide library file.");
    }
     double energy;
     while(nuclideFile >> energy){
        counts.push_back(energy);
     }
     nuclideFile.close();

    NuclidLibrary nuclideLibrary;
     nuclideLibrary.loadLibrary(nuclide_lib_path);
    std::vector<Nuclid> expected_nuclides = nuclideLibrary.getNuclideList();

    std::vector<double> baseline = BinBaseline::estimate(counts);

    this->ctx = BibSearchCtx{expected_nuclides, fwhm_model, ec_model, counts, baseline, Q0_conv_tol, iter};

}

std::vector<Nuclid> GammaMBibSearchDriver::run(const double chi2_threshold){


    std::vector<std::pair<Nuclid, FitOut>> results = GammaMBibSearch::search(ctx);

    std::vector<Nuclid> detected_nuclides;

    for (const auto& [nucl, fit] : results) {
        if (fit.reduced_chi2 <= chi2_threshold) {
            detected_nuclides.push_back(nucl);
        }
    }

    return detected_nuclides;
}