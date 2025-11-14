//
// Created by f.willems on 20.10.2025.
//

#pragma once

struct Bin {
    float bin_width;
    float min_x;
    float center_x;
    double count_average;

    Bin(double min_x, double bin_width){
        this->min_x = min_x;
        this->bin_width = bin_width;
        this->center_x = min_x + bin_width * 0.5;
        this->count_average = 0.0;
    }
};
