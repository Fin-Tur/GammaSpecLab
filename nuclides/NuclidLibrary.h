#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "Nuclid.h"

class NuclidLibrary {

    public:
    NuclidLibrary();
    void loadLibrary(const std::string& filepath);
    std::vector<Nuclid> getNuclideList() const;
    std::vector<Nuclid> getNuclidesNearEnergy(double energy_keV, double tolerance_keV) const;
    std::optional<Nuclid> getNuclide(int atomic_number, int mass_number) const;
    void buildNuclideMap();

    private:
    std::unordered_map<std::pair<int, int>, Nuclid> nuclide_map;

};