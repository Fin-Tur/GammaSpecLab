#include "API.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "../baseline/BinBaseline.h"

int estimate_baseline_bins(const double* counts, std::size_t count_size, double* baseline) {
    if (!counts || !baseline || count_size == 0) return -1;

    try {
        std::vector<double> input(counts, counts + count_size);
        std::vector<double> estimate = BinBaseline::estimate(input);
        if (estimate.size() != count_size) return -2;

        std::copy(estimate.begin(), estimate.end(), baseline);
        return 0;
    } catch (const std::exception&) {
        return -3;
    } catch (...) {
        return -4;
    }
}
