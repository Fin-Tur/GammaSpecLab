#ifndef TOOLBOX_EC_H
#define TOOLBOX_EC_H

struct EC {
    double offset;
    double slope;
    double quad;

    inline double energy_at(int channel) const {
        return offset + channel * slope + quad * channel * channel;
    }

    inline double dEdk(int channel) const{
        double gradient_k = slope + 2.0 * quad * channel;
        return gradient_k;
    }
};

#endif //TOOLBOX_EC_H