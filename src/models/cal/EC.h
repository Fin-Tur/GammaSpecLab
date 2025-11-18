#include <cmath>
#ifndef TOOLBOX_EC_H
#define TOOLBOX_EC_H

struct EC {
    double offset;
    double slope;
    double quad;

    //EC(double offset, double slope, double quad) : offset(offset), slope(slope), quad(quad) {}

    double energy_at(int channel) const {
        return offset + channel * slope + quad * channel * channel;
    }

    double dEdk(int channel) const{
        double gradient_k = slope + 2.0 * quad * channel;
        return gradient_k;
    }

    //EC only quadratic, so quadratic formula finds the channel- if deg(EC)>2, use binary search or approximation
    int channel_at(double energy) const {

        if(std::abs(quad) < 1e-12){
            if(std::abs(slope) < 1e-12){
                return static_cast<int>(offset + energy); //constant
            }
            return static_cast<int>((energy - offset) / slope);//linear
        }

        double a = quad;
        double b = slope;
        double c = offset - energy;
        
        double k1 = (-b + std::sqrt(b*b -4*a*c)) / (2.0 * a);
        double k2 = (-b - std::sqrt(b*b -4*a*c)) / (2.0 * a);

        if(k1 >=0 && k2 >= 0){
            return static_cast<int>(std::min(k1, k2));
        }if(k1 >=0){
            return static_cast<int>(k1);
        }if(k2 >=0){
            return static_cast<int>(k2);
        }
        //error vals
        return -1;
    }
};

#endif //TOOLBOX_EC_H