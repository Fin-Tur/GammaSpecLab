#ifndef TOOLBOX_FWHMC_H
#define TOOLBOX_FWHMC_H

struct FWHMC
{
    double offset;
    double slope;
    double quad;

    FWHMC(double offset, double slope, double quad) : offset(offset), slope(slope), quad(quad) {}

    inline double fwhm_at(double energy) const {
        return offset + energy * slope + quad * energy * energy;
    }

    inline double fwhm_at_channel(int channel) const{
        double return_dummy = 1; //FIXME imlementation fwhm channel, or conversion from fwhm_at (energy)
        return return_dummy;
    }

};


#endif //TOOLBOX_FWHMC_H