#ifndef TOOLBOX_GAMMABIBSEARCH_H
#define TOOLBOX_GAMMABIBSEARCH_H

#include <vector>
#include "../models/elements/FitOut.h"
#include "../models/ctx/BibSearchCtx.h"

namespace GammaMBibSearch {

    double calculate_spectral_gradient(const BibSearchCtx &ctx, int channel);
    double gaussian_base(const BibSearchCtx& ctx,int channel, int peak_center, double z);
    FitOut fit_once_energy(const BibSearchCtx &ctx, const std::vector<double>& Ycorr, int peak_center, double z_keV, double Ec_keV);
    FitOut fit_with_Q0_energy(const BibSearchCtx &ctx, const std::vector<double>& Ycorr, int peak_center, double z_keV);

    std::vector<std::pair<Nuclid, FitOut>> search(BibSearchCtx& ctx);
};

#endif //TOOLBOX_GAMMABIBSEARCH_H