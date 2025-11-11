#ifndef TOOLBOX_GAMMABIBSEARCH_H
#define TOOLBOX_GAMMABIBSEARCH_H

#include <vector>
#include "../models/Bin.h"
#include "../models/FitOut.h"
#include "../models/LinearBGCtx.h"
#include "../models/EC.h"
#include "../models/BibSearchCtx.h"

namespace GammaMBibSearch {

    LinearBGCtx make_linear_bg_ctx(const BibSearchCtx &ctx, int lR, int rR, int n);
    double calc_y_linear(const BibSearchCtx& ctx, int channel);
    double sigma2_genie_linear_bg(const BibSearchCtx& ctx, double Ycorr, int i);
    double calculate_spectral_gradient(const BibSearchCtx &ctx, int channel);
    double gaussian_base(const BibSearchCtx& ctx,int channel, int peak_center, double z);
    double gaussian_simple(const BibSearchCtx& ctx,int channel, int peak_center, int counts_center, double z);
    double gaussian_corrected(const BibSearchCtx& ctx, int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0);
    double gaussian_corrected_with_offset(const BibSearchCtx& ctx, int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double B0, double B1);

    FitOut fit_once_energy(const BibSearchCtx& ctx, int peak_center, double z_keV, double Ec_keV);
    FitOut fit_with_Q0_energy(const BibSearchCtx& ctx, int peak_center, double z_keV);

    std::vector<FitOut> search(const BibSearchCtx& ctx);
};

#endif //TOOLBOX_GAMMABIBSEARCH_H