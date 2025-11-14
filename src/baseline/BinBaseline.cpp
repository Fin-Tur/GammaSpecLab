//
// Created by f.willems on 13.10.2025.
//

//TODO: Adding parameters for FWHM; rn in var_needed.h

#include "BinBaseline.h"

#include <cmath>
#include "../models/FWHMC.h"

static thread_local FWHMC fwhm_model;

std::vector<double> BinBaseline::estimate(const std::vector<double> &counts) {

    //Magic Numbers -> Genie Manual :
    float k = 0.4;  
    int max_channels = 3400;
    int iters = 8;
    int number_of_plys = 4;
    //====================================================

    int n_counts = counts.size();
    if(n_counts <= 0) return counts;
    auto f1   = std::max(1., (double)fwhm_model.fwhm_at(1));
    auto fmid = std::max(1., (double)fwhm_model.fwhm_at(std::max(1, n_counts/2)));

    int number_of_bins = n_counts * f1 / (fmid * k);

    while (number_of_bins > max_channels) {
        k *= 1.1;
        number_of_bins = n_counts * f1 / ((fmid) * k);
    }

    //Create & Adjust bins & PLYs
    std::vector<std::vector<Bin>> plys(number_of_plys);
    std::vector<Bin> bins(number_of_bins);
    bins[0].min_x = 0;
    bins[0].bin_width = k * fwhm_model.fwhm_at(1);
    bins[0].center_x = bins[0].min_x + bins[0].bin_width * 0.5;
    plys[0].emplace_back(bins[0]);

    for (int i = 1; i < number_of_bins; i++) {
        Bin& bin = bins[i];
        bin.min_x = bins[i-1].min_x + bins[i-1].bin_width;
        bin.bin_width = std::max(k * fwhm_model.fwhm_at(bin.min_x), 1.);
        bin.center_x = bin.min_x + bin.bin_width * 0.5;
        //Calculate median channel counts
        double sum_counts = 0;
        double n = 0;
        for (int j = std::floor(bin.min_x); j < std::ceil(bin.min_x+bin.bin_width); ++j) {
            if (j >= n_counts) break;
            if (j - bin.min_x < 0) {
                double w = 1 -(bin.min_x - j);
                sum_counts += counts[j] * w;
                n += w;
            }else if (j + 1.0 > bin.min_x+bin.bin_width) {
                double w = (bin.min_x+bin.bin_width) - j;
                sum_counts += counts[j] * w;
                n += w;
            }else {
                sum_counts += counts[j];
                n+=1;
            }

        }
        if (n > 0) {
            bin.count_average = sum_counts / n;
        } else {
            bin.count_average = 0;
        }
        plys[i%number_of_plys].emplace_back(bins[i]);
    }

    //Estimating Baseline
    //------------------------------------------------
    //Since number_of_plys = 4 : 
    //Genie Manusal does not describe edge cases, since B'k is = Bk v (B'k-4 + B'k+4) /2
    //We simply skip first and last bin of each ply ?/ not a good solution but its not mentioned anywhere
    //-------------------------------------------------
    for(int i = 0; i < iters; i++) {
        for(auto& ply : plys){
            for(int i = 1; i < ply.size()-1; i++) {
                Bin& bin = ply[i];
                double prev_baseline = ply[i-1].count_average;
                double next_baseline = ply[i+1].count_average;
                bin.count_average = std::min(bin.count_average, (prev_baseline + next_baseline) * 0.5);
            }
        }
    }


    //Map baselines back to counts using linear interpolation
    //Edge cases not handled yet
    std::vector<double> counts_baseline(n_counts, 0.0);
    double E_i;
    for(int i = 0; i < n_counts; i++) {
        E_i = std::numeric_limits<double>::infinity();
        for(auto& ply : plys) {
            Bin* bin_left = nullptr;
            Bin* bin_right = nullptr;
            for(int j = 0; j < ply.size()-1; j++) {
                if(i >= ply[j].center_x && i < ply[j+1].center_x) {
                    bin_left = &ply[j];
                    bin_right = &ply[j+1];
                    break;
                }
            }
            if(bin_left == nullptr || bin_right == nullptr) continue;
            if(bin_right->center_x - bin_left->center_x < 1e-6) continue;
            double slope = (bin_right->count_average - bin_left->count_average) / (bin_right->center_x - bin_left->center_x);
            E_i = std::min(E_i, bin_left->count_average + slope * (i - bin_left->center_x));
            delete bin_left;
            delete bin_right;
        }

        counts_baseline[i] = std::min(E_i, counts[i]);
    }


    return counts_baseline;

    //Continuum is respectively estimated lower than it actually is : Has to be corrected later on (East, L.V., Phillips, R.L. and Strong, A.R. (1982). Nucl. Instr. & Meth. 199)
}
