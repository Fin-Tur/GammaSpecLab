#include "./NuclidLibrary.h"
#include <fstream>
#include <filesystem>
#include <optional>

//File has to look like this: atomic_number mass_number half_life gamma_energy gamma_intensity \n ...
//Reading single attributes separated by spaces/backshlashes/tabs

void NuclidLibrary::loadLibrary(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Failed to open file");
    }
    
    Nuclid nuclide;
    while (file >> nuclide.atomic_number >> nuclide.mass_number >> nuclide.half_life >> nuclide.gamma_energy >> nuclide.gamma_intensity) {
        nuclide_map[{nuclide.atomic_number, nuclide.mass_number}] = nuclide;
    }
}

std::vector<Nuclid> NuclidLibrary::getNuclideList() const {
    return std::vector<Nuclid>(nuclide_map.begin(), nuclide_map.end());
}

std::vector<Nuclid> NuclidLibrary::getNuclidesNearEnergy(double energy_keV, double tolerance_keV) const {
    std::vector<Nuclid> res;
    for(const auto& [key, nuc] : this->nuclide_map){
        if(std::abs(nuc.gamma_energy - energy_keV) <= tolerance_keV){
            res.push_back(nuc);
        }
    }
    return res;
}

std::optional<Nuclid> NuclidLibrary::getNuclide(int atomic_number, int mass_number) const {
    if (auto it = nuclide_map.find({atomic_number, mass_number}); it != nuclide_map.end()) {
        return it->second;
    }
    return std::nullopt;
}