
#pragma once
#include <stddef.h>

#if defined(_WIN32)
  #define EXPORT __declspec(dllexport)
#elif defined(__GNUC__)
  #define EXPORT __attribute__((visibility("default")))
#else
  #define EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif
    //0 OK, -1 invalid args, -2 size mismatch, -3 runtime error, -4 unknown
    EXPORT int baseline_estimate_bins(const double* counts, size_t count_size, int bins, double* baseline_out);
    EXPORT int baseline_estimate_als(const double* counts, size_t count_size, double lambda, double p, int max_iterations, double* baseline_out);

#ifdef __cplusplus
}
#endif
