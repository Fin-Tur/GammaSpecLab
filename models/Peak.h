#include <cmath>
#include "../var_needed.h"

#ifndef PEAK_H
#define PEAK_H

struct Peak{
    int id;
    int counts_center;
    int peak_center;
    double sigma;
    double z;
    Peak(int id, int counts_center, int peak_center, double sigma, double z) : id(id), counts_center(counts_center), peak_center(peak_center), sigma(sigma), z(z) {}
};

#endif //PEAK_H