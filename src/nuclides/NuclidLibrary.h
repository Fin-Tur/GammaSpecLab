#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "../models/elements/Nuclid.h"
#include <optional>

class NuclidLibrary {

    public:
    NuclidLibrary() = default;
    void loadLibrary(const std::string& filepath);
    std::vector<Nuclid> getNuclideList() const;
    std::vector<Nuclid> getNuclidesNearEnergy(double energy_keV, double tolerance_keV) const;
    std::optional<Nuclid> getNuclide(int atomic_number, int mass_number) const;

    private:
    std::vector<Nuclid> nuclides;

};