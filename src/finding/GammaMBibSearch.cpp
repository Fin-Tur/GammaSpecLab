#include "GammaMBibSearch.h"
#include <cmath>
#include "../thirdparty/eigen/Eigen/Sparse"
#include "../thirdparty/eigen/Eigen/QR"
#include "../baseline/BinBaseline.h"
#include "../models/cal/FWHMC.h"
#include "../models/ctx/BibSearchCtx.h"
using Eigen::MatrixXd; using Eigen::VectorXd;

//Calculates spectral gradient of spectrum at given channel
double GammaMBibSearch::calculate_spectral_gradient(const BibSearchCtx &ctx, int channel) {
    if(channel <= 2 || channel >= ctx.counts.size() - 2) return 0;
    //Genie Manuel: 0.5 * ((counts[channel + 2] + (!)1*counts[channel + 1])*0.25 - (2*counts[channel - 1] + counts[channel - 2])*0.25) - probably wrong so corrected to:
    double gradient_k = 0.5 * ((ctx.counts[channel + 2] + 2*ctx.counts[channel + 1])*0.25 - (2*ctx.counts[channel - 1] + ctx.counts[channel - 2])*0.25);
    return ctx.ec_model.dEdk(channel) != 0.0 ?  gradient_k / ctx.ec_model.dEdk(channel) : 0.0;
}

    
//Gaussian Base Model for Matrix Construction
double GammaMBibSearch::gaussian_base(const BibSearchCtx& ctx,int channel, int peak_center, double z) {
    double exponent = -std::pow(ctx.ec_model.energy_at(channel) - ctx.ec_model.energy_at(peak_center), 2) / (z * z);
    return std::exp(exponent); 
}

//1 FitSchritt 1 Peak
FitOut GammaMBibSearch::fit_once_energy(const BibSearchCtx &ctx, const std::vector<double>& Ycorr, int peak_center, double z_keV, double Ec_keV){
    const int N = ctx.rR - ctx.lR +1;
    MatrixXd A(N,4);
    VectorXd y(N), w(N);

    for(int row = 0, k = ctx.lR; k <= ctx.rR; ++k, ++row){
       const double Ei = ctx.ec_model.energy_at(k);
       const double grad = calculate_spectral_gradient(ctx, k); 
       const double phi = gaussian_base(ctx, k, peak_center, z_keV);

       //Matrix Columns:
       A(row,0) = 1.0;                //b0
       A(row,1) = (Ei - Ec_keV);      //b1
       A(row,2) = grad;               //dY/dE
       A(row,3) = phi;                //counts_center*Gauss(E)

       y(row) = Ycorr[k];
       //const double sigma2 = sigma2_genie_linear_bg(ctx, Ycorr[k], k); //FIXME
       const double sigma2 = z_keV;
       w(row) = (sigma2 > 0.0) ? 1.0 / sigma2 : 0.0;
    }

    //Weighted Least Squares Solution: x = (A^T W A)theta A^T W y
    const VectorXd w_sqrt = w.array().sqrt();
    const auto w_diag = w_sqrt.asDiagonal();
    auto Ahat = w_diag * A;
    const auto yhat = w_diag * y;

    Eigen::ColPivHouseholderQR<MatrixXd> qr(Ahat);
    VectorXd theta = qr.solve(yhat);

    VectorXd residuals = y - A * theta;

    const double chi2 = residuals.transpose() * (w.asDiagonal() * residuals);
    int dof = A.rows()- A.cols();
    const double reduced_chi2 = (dof > 0) ? chi2 / dof : 1.0;

    return {
        theta(0), //b0
        theta(1), //b1
        theta(2), //Q0
        theta(3),  //counts_center
        reduced_chi2 //fit quality
    };
}

FitOut GammaMBibSearch::fit_with_Q0_energy(const BibSearchCtx &ctx, const std::vector<double>& Ycorr, int peak_center, double z_keV){
    double Ec = ctx.ec_model.energy_at(peak_center);
    double last_Q0 = 0.0;
    FitOut fit_result{};

    for(int it = 0; it < 8; ++it){
        fit_result = fit_once_energy(ctx, Ycorr, peak_center, z_keV, Ec);
        const double dQ0 = fit_result.Q0 - last_Q0;
        if(std::abs(dQ0) < ctx.Q0_conv_tol) break;
        Ec += fit_result.Q0;
        last_Q0 = fit_result.Q0;
    }
    fit_result.z_keV = z_keV;
    return fit_result;
}

std::vector<std::pair<Nuclid, FitOut>> GammaMBibSearch::search(BibSearchCtx& ctx) {

    std::vector<std::pair<Nuclid, FitOut>> detected_peaks;

    std::vector<double> corrected_spectrum(ctx.counts.size());
    std::transform(ctx.counts.begin(), ctx.counts.end(), ctx.baseline.begin(), corrected_spectrum.begin(), [](double count, double base) { return count - base; });

    for(auto& nucl : ctx.expected_nuclides) {

        const int ck = ctx.ec_model.channel_at(nucl.gamma_energy);

        const double Ec = ctx.ec_model.energy_at(ck);
        const double FWHM_keV = ctx.fwhm_model.fwhm_at(Ec);
        const double z_keV = FWHM_keV / (std::sqrt(2.0*std::log(2)));

        const int N = std::ceil((z_keV/ctx.ec_model.dEdk(ck))*1.5);
        const int lR = std::max(0, ck - N/2);
        const int rR = std::min(ck + N/2, static_cast<int>(ctx.counts.size() - 1));

        ctx.lR = lR;
        ctx.rR = rR;

        auto fit = fit_with_Q0_energy(ctx, corrected_spectrum, ck, z_keV);
        
        detected_peaks.emplace_back(nucl, fit);

    }


    return detected_peaks;
}
