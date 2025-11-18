#pragma once
#include <vector>
#include "../cal/EC.h"
#include "../cal/FWHMC.h"
#include "../elements/Peak.h"

struct SecDervCtx {

    std::vector<double> counts;
    FWHMC fwhm_cal;
    EC ec_cal;

    float s_min; //Genie Manuel recommends > 3
    int blocksize;

    //to be filled
    std::vector<Peak> found_peaks;
};