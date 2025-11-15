# GammaDecayLib

A C++ library for scientific analysis of gamma-decay spectra. This
project provides modular components for baseline estimation, peak
finding, nuclide matching, and calibration workflows. It integrates
numerical linear algebra via Eigen and exposes a unified API for
spectrum-processing pipelines.

## Overview

The library is structured into independent modules:

### API

Future accespoint to various algorithmns, accessible via .dll build. Algorithms can be used independently, or ready-to-use Drivers.

### Baseline

Implements algorithms such as Asymmetric Least Squares (ALS) and binned
baseline estimation.

### Peak Finding

Includes second-derivative based search and a multi-bibliography
gamma-line search engine.

### Peak Fitting

Includes basic gaussian models as part of Gamma-M Library search, 
soon to be extended with advanced models using Taylor Series to linearize variables in the fitting algorithmn

### Models

Core data structures: - **cal**: Energy and shape calibration utilities .
                      - **ctx**: overall use of context based data structures to keep workflows organized.

### Nuclides

Nuclide database utilities for mapping detected peaks to isotopic
candidates. Nuclid databases are not given atp and have to be included by the user.

### Drivers

Executable drivers for testing peak search algorithms and bibliography
matching.

### Third‑Party

Bundled Eigen for linear algebra computations.

## Build

This project uses CMake:

``` bash
mkdir build && cd build
cmake ..
make
```

## Intended Use

The library is designed for:

-   Gamma spectrum preprocessing
-   Baseline estimation
-   Peak detection through derivative-informed and bibliographic methods
-   Mapping detected peaks to nuclides via curated decay tables
-   Integration into larger spectroscopy toolchains

## Future Extensions

Planned extensions include: - Advanced calibration models -Full-spectrum fitting
and uncertainty propagation - API via .dll / .so
