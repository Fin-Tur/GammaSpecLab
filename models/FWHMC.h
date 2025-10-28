#ifndef TOOLBOX_FWHMC_H
#define TOOLBOX_FWHMC_H

struct FWHMC
{
    double offset;
    double slope;
    double quad;

    inline double fwhm_at(double energy) const {
        return offset + energy * slope + quad * energy * energy;
    }

};


#endif //TOOLBOX_FWHMC_H