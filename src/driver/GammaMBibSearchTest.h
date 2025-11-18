//
// Created by f.willems on 18.11.2025.
//

#ifndef TOOLBOX_GAMMAMBIBSEARCHTEST_H
#define TOOLBOX_GAMMAMBIBSEARCHTEST_H
#include "../models/ctx/BibSearchCtx.h"

class GammaMBibSearchTest{
public:
    GammaMBibSearchTest();
    bool test_gamma_m_bib_search();
private:
    BibSearchCtx ctx;
};

#endif //TOOLBOX_GAMMAMBIBSEARCHTEST_H