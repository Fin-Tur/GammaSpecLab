#ifndef TOOLBOX_GAMMABIBSEARCH_H
#define TOOLBOX_GAMMABIBSEARCH_H

#include <vector>
#include "../var_needed.h"
#include "../models/Bin.h"
#include "../models/Peak.h"
#include "../models/FitOut.h"
#include "../models/LinearBGCtx.h"
#include "../models/EC.h"

namespace GammaMBibSearch {

    LinearBGCtx make_linear_bg_ctx(const std::vector<double> &counts, int lR, int rR, int n);
    double calc_y_linear(const std::vector<double> &counts, LinearBGCtx &bg_ctx, int channel);
    double sigma2_genie_linear_bg(int i);
    double calculate_spectral_gradient(const std::vector<double> &counts, int channel);
    double gaussian_base(int channel, int peak_center, double z);
    double gaussian_simple(int channel, int peak_center, int counts_center, double z);
    double gaussian_corrected(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0);
    double gaussian_corrected_with_offset(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double b0, double b1);

    FitOut fit_once_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, double Ec_keV);
    FitOut fit_with_Q0_energy(const std::vector<double>& Ycorr, int lR, int rR, int peak_center, double z_keV, int iter = 8, double conv_tol = 1e-3);

    std::vector<FitOut> search(const std::vector<double> &counts);
};

#endif //TOOLBOX_GAMMABIBSEARCH_H