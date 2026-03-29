#include <fstream>
#include <string>

bool fichiersEquivalents(const std::string& fichier1, const std::string& fichier2) {
    std::ifstream f1(fichier1);
    std::ifstream f2(fichier2);

    // Vérifier que les fichiers sont bien ouverts
    if (!f1.is_open() || !f2.is_open()) {
        return false;
    }

    std::string ligne1, ligne2;

    while (true) {
        bool lire1 = static_cast<bool>(std::getline(f1, ligne1));
        bool lire2 = static_cast<bool>(std::getline(f2, ligne2));

        // Si les deux fichiers sont terminés en même temps
        if (!lire1 && !lire2) {
            return true;
        }

        // Si un fichier se termine avant l'autre ou lignes différentes
        if (lire1 != lire2 || ligne1 != ligne2) {
            return false;
        }
    }
}