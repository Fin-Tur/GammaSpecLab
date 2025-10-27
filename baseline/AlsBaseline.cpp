//
// Created by f.willems on 07.08.2025.
//
#include "thirdparty/eigen/Eigen/Sparse"

#include "AlsBaseline.h"

//Estimate Background using ALS
std::vector<double> ALS_BASELINE::estimateBackgroundUsingALS(std::vector<double> counts, double lambda, double p, int maxIterations) {
    //initialize len, weights
    int cntLen = counts.size();
    std::vector<double> weights(cntLen, 1.0);

    //Create second derivation Matrix / curvature
    Eigen::SparseMatrix<double> D(cntLen - 2, cntLen);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < cntLen - 2; ++i) {
        triplets.emplace_back(i, i, 1.0);
        triplets.emplace_back(i, i + 1, -2.0);
        triplets.emplace_back(i, i + 2, 1.0);
    }
    D.setFromTriplets(triplets.begin(), triplets.end());

    //Initialize Curvature Penalty Matrix
    Eigen::SparseMatrix<double> CurvaturePenalty = lambda * (D.transpose() * D);
    //Initialize counts vector and background vec
    Eigen::VectorXd counts_Vec = Eigen::Map<const Eigen::VectorXd>(counts.data(), counts.size());
    Eigen::VectorXd background = counts_Vec;
    //Fit background noise for iter iterations
    for (int iter = 0; iter < maxIterations; ++iter) {
        //Create and fill weighted matrix
        std::vector<Eigen::Triplet<double>> Weighted_Diagonal_Sparse_Matrix;
        for (int i = 0; i < cntLen; ++i) {
            Weighted_Diagonal_Sparse_Matrix.emplace_back(i, i, weights[i]);
        }
        Eigen::SparseMatrix<double> Weights(cntLen, cntLen);
        Weights.setFromTriplets(Weighted_Diagonal_Sparse_Matrix.begin(), Weighted_Diagonal_Sparse_Matrix.end());
        // A = weight + penalty
        Eigen::SparseMatrix<double> A = Weights + CurvaturePenalty;
        //Calculates weight * count
        Eigen::VectorXd Weighted_counts = counts_Vec.cwiseProduct(Eigen::Map<Eigen::VectorXd>(weights.data(), cntLen));

        //Solve LGS to estimate backgrounds
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        background = solver.solve(Weighted_counts);

        //fit weights and break off if changed weight is lower then 1e-6
        double delta = 0.0;
        for (int i = 0; i < cntLen; ++i) {
            double diff = counts[i] - background[i];
            double w_before = weights[i];
            weights[i] = diff > 0 ? p : 1.0 - p;
            double w_after = weights[i];
            delta += std::abs(w_before - w_after);
        }

        if (delta < 1e-6) break;
    }

    //Return estimated Background
    return std::vector<double>(background.data(), background.data() + background.size());
}

