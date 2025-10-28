#include "GammaMBibSearch.h"
#include <cmath>
#include "thirdparty/eigen/Eigen/Sparse"
#include "thirdparty/eigen/Eigen/QR"
#include "baseline/BinBaseline.h"
#include <FWHMC.h>


using Eigen::MatrixXd; using Eigen::VectorXd;

static thread_local LinearBGCtx last_bg_ctx;
//For further processing, we will expect this vector to be filled with our library peaks
static thread_local std::vector<Peak> expected_peaks;
static thread_local FWHMC fwhm_model{FWHM_OFFSET, FWHM_SLOPE, FWHM_QUAD};
static thread_local EC ec_model{EC_OFFSET, EC_SLOPE, EC_QUAD};

//Calculates spectral gradient of spectrum at given channel
double GammaMBibSearch::calculate_spectral_gradient(const std::vector<double> &counts, int channel) {
    if(channel <= 2 || channel >= counts.size() - 2) return 0;
    //Genie Manuel: 0.5 * ((counts[channel + 2] + (!)1*counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25) - probably wrong so corrected to:
    double gradient_k = 0.5 * ((counts[channel + 2] + 2*counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25);
    return ec_model.dEdk(channel) != 0.0 ?  gradient_k / ec_model.dEdk(channel) : 0.0;
}

//calcs linear_bg_ctx
LinearBGCtx GammaMBibSearch::make_linear_bg_ctx(const std::vector<double> &counts, int lR, int rR, int n){
    int N = rR - lR +1;
    double B0 = 0.0, B1 = 0.0;
    for(int i = std::max(0, lR - n); i <= lR -1; ++i){
        B0 += counts[i];
    }   
    for(int i = rR +1; i <= std::min(static_cast<int>(counts.size()-1), rR + n); ++i){
        B1 += counts[i];
    }
    return {lR, rR, n, N, B0, B1};
}

//calcs y_linear after bg removal GM:(6) 
//eigentlich unnötig warum habe ich das überhaupt gesachrieben ??
double GammaMBibSearch::calc_y_linear(const std::vector<double> &counts, int channel){
    if(!(channel >= last_bg_ctx.lR && channel <= last_bg_ctx.rR)) {
        return counts[channel];
    }
    int i = channel-last_bg_ctx.lR + 1;
    double y_apostrophe = counts[channel] - (1.0/last_bg_ctx.n)*last_bg_ctx.B0 - (i*(last_bg_ctx.B1 - last_bg_ctx.B0))/(last_bg_ctx.n*(last_bg_ctx.N + 1));
    return y_apostrophe;
}


//sigma2 according to Genie Manual with linear bg GM:(7)
double GammaMBibSearch::sigma2_genie_linear_bg(const double Ycorr, int i) {
    int k = (i - last_bg_ctx.lR + 1);
    double topTerm = last_bg_ctx.B0*(std::pow(last_bg_ctx.N+1-k, 2)) + last_bg_ctx.B1*(k*k);
    double bottomTerm = std::pow(last_bg_ctx.n*(last_bg_ctx.N+1), 2);
    return bottomTerm != 0.0 ? Ycorr + (topTerm / bottomTerm) : Ycorr;
}

    
//Gaussian Base Model for Matrix Construction
double GammaMBibSearch::gaussian_base(int channel, int peak_center, double z) {
    double exponent = -std::pow(ec_model.energy_at(channel) - ec_model.energy_at(peak_center), 2) / (z * z);
    return std::exp(exponent); 
}

//simple Gaussian model without corrections
double GammaMBibSearch::gaussian_simple(int channel, int peak_center, int counts_center, double z) {
    double exponent = -std::pow(ec_model.energy_at(channel) - ec_model.energy_at(peak_center), 2) / (z * z);
    return (counts_center * std::exp(exponent));
}

//Gaussian corrected for spectral gradient
double GammaMBibSearch::gaussian_corrected(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0) {
    double exponent = -std::pow(ec_model.energy_at(channel) - ec_model.energy_at(peak_center), 2) / (z * z);
    double correction = (1 + Q0) * spec_gradient;
    return (counts_center * std::exp(exponent) + correction);
}

//Gaussian corrected for spectral gradient and baseline offset
double GammaMBibSearch::gaussian_corrected_with_offset(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double b0, double b1) {
    double exponent = -std::pow(ec_model.energy_at(channel) - ec_model.energy_at(peak_center), 2) / (z * z);
    double correction = (1 + Q0) * spec_gradient;
    double baseline = b0 + b1 * (ec_model.energy_at(channel) - ec_model.energy_at(peak_center));
    return (counts_center * std::exp(exponent) + correction + baseline);
}

//1 FitSchritt 1 Peak
FitOut GammaMBibSearch::fit_once_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, double Ec_keV){
    const int N = rR - lR +1;
    MatrixXd A(N,4);
    VectorXd y(N), w(N);

    for(int row = 0, k = lR; k <= rR; ++k, ++row){
       const double Ei = ec_model.energy_at(k);
       const double grad = calculate_spectral_gradient(Ycorr, k); 
       const double phi = gaussian_base(k, peak_center, z_keV);

       //Matrix Columns:
       A(row,0) = 1.0;                //b0
       A(row,1) = (Ei - Ec_keV);      //b1
       A(row,2) = grad;               //dY/dE
       A(row,3) = phi;                //counts_center*Gauss(E)

       y(row) = Ycorr[k];
       const double sigma2 = sigma2_genie_linear_bg(Ycorr[k], k); //TODO
       w(row) = (sigma2 > 0.0) ? 1.0 / sigma2 : 0.0;
    }

    //Weighted Least Squares Solution: x = (A^T W A)theta A^T W y
    const VectorXd w_sqrt = w.array().sqrt();
    MatrixXd Ahat = w_sqrt.asDiagonal() * A;
    VectorXd yhat = w_sqrt.asDiagonal() * y;

    Eigen::ColPivHouseholderQR<MatrixXd> qr(Ahat);
    VectorXd theta = qr.solve(yhat);

    double Q0 = theta(2) - 1.0; //dY/dE coefficient
    return {
        theta(0), //b0
        theta(1), //b1
        Q0,
        theta(3)  //counts_center
    };
}

FitOut GammaMBibSearch::fit_with_Q0_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, int iter, double conv_tol){
    double Ec = ec_model.energy_at(peak_center);
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

    //Parameters: E(channel_to_energy(peak_center)), Z(peak_width*sqrt(2)) are given by the Library and Calibration
    //Fitting can be continued linearly from there

    std::vector<double> baseline = BinBaseline::estimate(counts);
    std::vector<double> corrected_spectrum(counts.size());
    std::transform(counts.begin(), counts.end(), baseline.begin(), corrected_spectrum.begin(), [](double count, double base) { return count - base; });

    for(auto& peak : expected_peaks) {

        int ck = peak.peak_center;
        double Ec = ec_model.energy_at(ck);
        double FWHM_keV = fwhm_model.fwhm_at(Ec);
        double z_keV = FWHM_keV / (std::sqrt(2.0*std::log(2)));

        int N = std::ceil((z_keV/ec_model.dEdk(ck))*1.5);
        int lR = std::max(0, ck - N/2);
        int rR = std::min(ck + N/2, static_cast<int>(counts.size() - 1));

        //TODO: Dont know proper n yet
        int n = 1;
        //Important: lR-n and rR+n must be valid indices

        last_bg_ctx = make_linear_bg_ctx(corrected_spectrum, lR, rR, n);

        auto fit = fit_with_Q0_energy(corrected_spectrum, lR, rR, ck, z_keV);
        detected_peaks.emplace_back(fit);

    }


    return detected_peaks;
}
