#include <vector>
#include "LinearBGCtx.h"
#include "FWHMC.h"
#include "EC.h"
#include "FitOut.h"
#include "Nuclid.h"

#pragma once

struct BibSearchCtx{

    //set
    std::vector<Nuclid> expected_nuclides;
    FWHMC fwhm_model;
    EC ec_model;
    std::vector<double> counts;

    //user dependent
    double Q0_conv_tol;
    int iter;

    //empty on initialize- filled within search
    LinearBGCtx bg_ctx;
    //FitOut last_fit;
};