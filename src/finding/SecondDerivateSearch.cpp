#include "SecondDerivateSearch.h"

#include <numeric>
#include <cmath>

std::vector<Peak> SecondDerivateSearch::find_peaks(const SecDervCtx& ctx) {

    std::vector<Peak> peaks;

    //calc every (blocksize) channels
    int n_mid;
    double FWHM;
    double cw;
    int k;
    std::vector<double> coeff;

    //clac every channel
    double dd;
    double sd;
    double ss;

    std::vector<std::pair<int, double>> ss_values; //channel, ss

    for(int channel = 0; channel < ctx.counts.size(); ++channel){

        //every (blocksize) channels recalc cw and coeff
        if(channel % ctx.blocksize == 0){
            //calc cw
            n_mid = channel + ctx.blocksize /2;
            FWHM = ctx.fwhm_cal.fwhm_at_channel(std::min(n_mid, static_cast<int>(ctx.counts.size()-1)));
            cw = FWHM / 2.355;
            //calc coeff, k
            coeff.clear();
            k = 0;
            coeff.emplace_back(-100);

            for(int j = 1; std::abs(coeff[j-1]) >= 1; ++j){ 
                double c = (100*(j*j - (cw*cw)) / (cw*cw)) * std::exp(-(j*j)/(2*cw*cw));
                coeff.emplace_back(c);
                k = j;
            }
            coeff.resize(coeff.size()-1);
            k--;
            //normalize coeff[1] so sum(coeff)=0
            double sum = coeff[0] + 2.0 * std::accumulate(coeff.begin() +1, coeff.end(), 0.0);
            coeff[1] -= sum/2;
        }

        //calc dd and sd for every channel -> //TODO what to do for channel < k || counts_size -1 -k < channel 
        if(channel < k || ctx.counts.size() -1 -k < channel) continue;

        dd = 0.0;
        sd = 0.0;

        for(int j = -k; j <= k; ++j){
            double coeff_j = coeff[std::abs(j)]; 
            dd += coeff_j * ctx.counts[channel + j];
            sd += coeff_j * coeff_j * ctx.counts[channel + j];
            
        }
        ss = dd / std::sqrt(sd);
        if(ss < 0 && std::abs(ss) >= ctx.s_min){
            ss_values.emplace_back(channel, ss);
        }else{
            if(!ss_values.empty()){
                //construct peak
                peaks.emplace_back(construct_peak(ctx, ss_values));
                ss_values.clear();
            }
        }        

    }

    if(!ss_values.empty()){
        peaks.emplace_back(construct_peak(ctx, ss_values));
        ss_values.clear();
    }
    return peaks;
}

    Peak SecondDerivateSearch::construct_peak(const SecDervCtx& ctx , const std::vector<std::pair<int, double>>& ss_vals) {
        double c;
        double sum_top = 0.0;
        double sum_btm = 0.0;
        for(const auto& [channel, ss] : ss_vals) {
            sum_top += channel * ss;
            sum_btm += ss;
        }
        c = (sum_btm != 0.0) ? sum_top / sum_btm : -1.0;
        const double peak_height = ctx.counts[static_cast<int>(c)];
        const double FWHM_keV = ctx.fwhm_cal.fwhm_at(ctx.ec_cal.energy_at(c));
        const double z_keV = FWHM_keV / (std::sqrt(2.0*std::log(2)));
        return Peak{c, peak_height, z_keV};
    }