
#ifndef ALS_BASELINE_H
#define ALS_BASELINE_H

#include <vector>

namespace ALS_BASELINE {
    std::vector<double> estimateBackgroundUsingALS(std::vector<double> counts, double lambda, double p, int maxIterations);
}

#endif //ALS_BASELINE_H