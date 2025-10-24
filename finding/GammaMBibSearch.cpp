#include "GammaMBibSearch.h"
#include <cmath>
#include "thirdparty/eigen/Eigen/Sparse"
#include "baseline/BinBaseline.h"

//Calculates spectral gradient of spectrum at given channel
double GammaMBibSearch::calculate_spectral_gradient(const std::vector<double> &counts, int channel) {
    if(channel <= 2 || channel >= counts.size() - 2) return 0;
    return 0.5 * ((counts[channel + 2] + counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25);
}
    

//simple Gaussian model without corrections
double GammaMBibSearch::gaussian_simple(int channel, int peak_center, int counts_center, double z) {
    double exponent = -std::pow(TOOLBOX_VAR_NEEDED_H::channel_to_energy(channel) - TOOLBOX_VAR_NEEDED_H::channel_to_energy(peak_center), 2) / (2 * z * z);
    return (counts_center * std::exp(exponent));
}

//Gaussian corrected for spectral gradient
double GammaMBibSearch::gaussian_corrected(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0) {
    double exponent = -std::pow(TOOLBOX_VAR_NEEDED_H::channel_to_energy(channel) - TOOLBOX_VAR_NEEDED_H::channel_to_energy(peak_center), 2) / (2 * z * z);
    double correction = (1 + Q0) * spec_gradient;
    return (counts_center * std::exp(exponent) + correction);
}

//Gaussian corrected for spectral gradient and baseline offset
double GammaMBibSearch::gaussian_corrected_with_offset(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double b0, double b1) {
    double exponent = -std::pow(TOOLBOX_VAR_NEEDED_H::channel_to_energy(channel) - TOOLBOX_VAR_NEEDED_H::channel_to_energy(peak_center), 2) / (2 * z * z);
    double correction = (1 + Q0) * spec_gradient + b0 + b1 * channel;
    double baseline = b0 + b1 * (TOOLBOX_VAR_NEEDED_H::channel_to_energy(channel) - TOOLBOX_VAR_NEEDED_H::channel_to_energy(peak_center));
    return (counts_center * std::exp(exponent) + correction + baseline);
}

std::vector<Peak> GammaMBibSearch::search(const std::vector<double> &counts) {
    std::vector<Peak> detected_peaks;
    //For further processing, we will expect this vector to be filled with our library peaks
    std::vector<Peak> expected_peaks;

    // Parameters: E(channel_to_energy(peak_center)), Z(peak_width*sqrt(2)) are given by the Library and Calibration
    // Fitting can be continued linearly from there

    std::vector<double> baseline = BinBaseline::estimate(counts);
    std::vector<double> corrected_spectrum(counts.size());
    std::transform(counts.begin(), counts.end(), baseline.begin(), corrected_spectrum.begin(), [](double count, double base) { return count - base; });

    for(auto& peak : expected_peaks) {
        int lR = -1; //TODO: determine left range based on peak properties
        int rR = 1;  //TODO: determine right range based on peak properties

        for(int i = lR; i < rR; i++) {
            //Standart values
            int b0 = 1;
            int b1 = i;
            int Q0 = calculate_spectral_gradient(corrected_spectrum, i);
            double counts_center = counts[peak.peak_center];

            double observed_value = corrected_spectrum[peak.peak_center + i];
            std::vector<double> weights;
            for(int j = lR; j < rR; j++) {
                double weight = TOOLBOX_VAR_NEEDED_H::FWHM(j)/2*std::sqrt(2*std::log(2));
                weights.push_back(weight);
            }

            //Building Designer Matrix

    }

    


    return detected_peaks;
}
