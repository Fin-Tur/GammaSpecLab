#include <vector>
#include "LinearBGCtx.h"
#include "../cal/FWHMC.h"
#include "../cal/EC.h"
#include "../elements/FitOut.h"
#include "../elements/Nuclid.h"

#pragma once

struct BibSearchCtx{
    public:
    //BibSearchCtx();

    //set
    std::vector<Nuclid> expected_nuclides;
    FWHMC fwhm_model;
    EC ec_model;
    std::vector<double> counts;
    std::vector<double> baseline;

    //user dependent
    double Q0_conv_tol;
    int iter;

    //empty on initialize- filled within search renewed each iter
    int lR;
    int rR;
};