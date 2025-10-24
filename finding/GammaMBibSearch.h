#ifndef TOOLBOX_GAMMABIBSEARCH_H
#define TOOLBOX_GAMMABIBSEARCH_H

#include <vector>
#include "../var_needed.h"
#include "../models/Bin.h"
#include "../models/Peak.h"

namespace GammaMBibSearch {
    double calculate_spectral_gradient(const std::vector<double> &counts, int channel);
    double gaussian_simple(int channel, int peak_center, int counts_center, double z);
    double gaussian_corrected(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0);
    double gaussian_corrected_with_offset(int channel, int peak_center, int counts_center, double z, double spec_gradient, double Q0, double b0, double b1);

    std::vector<Peak> search(const std::vector<double> &counts);
};

#endif //TOOLBOX_GAMMABIBSEARCH_H