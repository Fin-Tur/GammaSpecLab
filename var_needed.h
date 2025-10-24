//
// Created by f.willems on 20.10.2025.
//

#ifndef TOOLBOX_VAR_NEEDED_H
#define TOOLBOX_VAR_NEEDED_H

constexpr double FWHM_OFFSET = 5;
constexpr double FWHM_SLOPE = 15;
constexpr double FWHM_QUAD = 30;

constexpr double channel_to_energy(const int channel) {
    return FWHM_OFFSET + channel * FWHM_SLOPE + FWHM_QUAD * channel * channel;
}

constexpr double FWHM(const int channel) {
    //TODO
    return 1;
}


#endif //TOOLBOX_VAR_NEEDED_H