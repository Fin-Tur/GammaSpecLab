#ifndef TOOLBOX_FITOUT_H
#define TOOLBOX_FITOUT_H

    struct FitOut { 

        double b0; 
        double b1; 
        double Q0; 
        double counts_center; //FIXME what?
        double z_keV;

        double reduced_chi2;


    };

#endif //TOOLBOX_FITOUT_H