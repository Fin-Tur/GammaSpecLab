#include "GammaMBibSearch.h"
#include <cmath>
#include "thirdparty/eigen/Eigen/Sparse"
#include "thirdparty/eigen/Eigen/QR"
#include "baseline/BinBaseline.h"

using Eigen::MatrixXd; using Eigen::VectorXd;

//Calculates spectral gradient of spectrum at given channel
double GammaMBibSearch::calculate_spectral_gradient(const std::vector<double> &counts, int channel) {
    if(channel <= 2 || channel >= counts.size() - 2) return 0;
    //Genie Manuel: 0.5 * ((counts[channel + 2] + (!)1*counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25) - probably wrong so corrected to:
    double gradient_k = 0.5 * ((counts[channel + 2] + 2*counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25);
    return TOOLBOX_VAR_NEEDED_H::dEdk(channel) != 0.0 ?  gradient_k / TOOLBOX_VAR_NEEDED_H::dEdk(channel) : 0.0;
}

//Genie-Varianz after linear BG-Removal -> σ_i^2  (TODO: (6)–(7) Genie Manual)
double GammaMBibSearch::sigma2_genie_linear_bg(int /*i*/) { return 1.0; }  // => w_i = 1/σ_i^2

    
//Gaussian Base Model for Matrix Construction
double GammaMBibSearch::gaussian_base(int channel, int peak_center, double z) {
    double exponent = -std::pow(energy_at(channel) - energy_at(peak_center), 2) / (z * z);
    return std::exp(exponent); 
}

//simple Gaussian model without corrections
double GammaMBibSearch::gaussian_simple(int channel, int peak_center, int counts_center, double z) {
    double exponent = -std::pow(energy_at(channel) - energy_at(peak_center), 2) / (z * z);
    return (counts_center * std::exp(exponent));
}

//Gaussian corrected for spectral gradient
double GammaMBibSearch::gaussian_corrected(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0) {
    double exponent = -std::pow(energy_at(channel) - energy_at(peak_center), 2) / (z * z);
    double correction = (1 + Q0) * spec_gradient;
    return (counts_center * std::exp(exponent) + correction);
}

//Gaussian corrected for spectral gradient and baseline offset
double GammaMBibSearch::gaussian_corrected_with_offset(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double b0, double b1) {
    double exponent = -std::pow(energy_at(channel) - energy_at(peak_center), 2) / (z * z);
    double correction = (1 + Q0) * spec_gradient + b0 + b1 * channel;
    double baseline = b0 + b1 * (energy_at(channel) - energy_at(peak_center));
    return (counts_center * std::exp(exponent) + correction + baseline);
}

//1 FitSchritt 1 Peak
FitOut GammaMBibSearch::fit_once_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, double Ec_keV){
    const int N = rR - lR +1;
    MatrixXd A(N,4);
    VectorXd y(N), w(N);

    for(int row = 0, k = lR; k <= rR; ++k, ++row){
       const double Ei = energy_at(k);
       const double grad = calculate_spectral_gradient(Ycorr, k); 
       const double phi = gaussian_base(k, peak_center, z_keV);

       //Matrix Columns:
       A(row,0) = 1.0;                //b0
       A(row,1) = (Ei - Ec_keV);      //b1
       A(row,2) = grad;               //Q0 * dY/dE
       A(row,3) = phi;                //counts_center*Gauss(E)

       y(row) = Ycorr[k];
       const double sigma2 = sigma2_genie_linear_bg(k); //TODO
       w(row) = (sigma2 > 0.0) ? 1.0 / sigma2 : 0.0;
    }

    //Weighted Least Squares Solution: x = (A^T W A)theta A^T W y
    const VectorXd w_sqrt = w.array().sqrt();
    MatrixXd Ahat = A * w_sqrt.asDiagonal();
    VectorXd yhat = y * w_sqrt.asDiagonal();

    Eigen::ColPivHouseholderQR<MatrixXd> qr(Ahat);
    VectorXd theta = qr.solve(yhat);

    return {
        theta(0), //b0
        theta(1), //b1
        theta(2), //Q0
        theta(3)  //counts_center
    };
}

FitOut GammaMBibSearch::fit_with_Q0_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, int iter, double conv_tol){
    double Ec = energy_at(peak_center);
    double last_Q0 = 0.0;
    FitOut fit_result{};

    for(int it = 0; it < 8; ++it){
        fit_result = fit_once_energy(Ycorr, lR, rR, peak_center, z_keV, Ec);
        const double dQ0 = fit_result.Q0 - last_Q0;
        if(std::abs(dQ0) < 1e-3) break; //Convergence test
        Ec += fit_result.Q0;
        last_Q0 = fit_result.Q0;
    }
    return fit_result;
}

std::vector<FitOut> GammaMBibSearch::search(const std::vector<double> &counts) {
    std::vector<FitOut> detected_peaks;
    //For further processing, we will expect this vector to be filled with our library peaks
    std::vector<Peak> expected_peaks;

    //Parameters: E(channel_to_energy(peak_center)), Z(peak_width*sqrt(2)) are given by the Library and Calibration
    //Fitting can be continued linearly from there

    std::vector<double> baseline = BinBaseline::estimate(counts);
    std::vector<double> corrected_spectrum(counts.size());
    std::transform(counts.begin(), counts.end(), baseline.begin(), corrected_spectrum.begin(), [](double count, double base) { return count - base; });

    for(auto& peak : expected_peaks) {
        int lR = -1; //TODO: determine left range based on peak properties
        int rR = 1;  //TODO: determine right range based on peak properties

        int ck = peak.peak_center;
        double Ec = energy_at(ck);
        double FWHM_keV = fwhm_at(Ec);
        double z_keV = FWHM_keV / (std::sqrt(2.0*std::log(2)));

        auto fit = fit_with_Q0_energy(corrected_spectrum, lR, rR, ck, z_keV);
        detected_peaks.emplace_back(fit);

    }


    return detected_peaks;
}
