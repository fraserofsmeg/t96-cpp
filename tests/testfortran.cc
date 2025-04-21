#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#define _USE_MATH_DEFINES

extern "C" {
    void t96_(int *iopt, double *parmod, double *ps,
               double *x, double *y, double *z,
               double *bx, double *by, double *bz);
}

void t96(int iopt, double* parmod, double ps,
         double x, double y, double z,
         double* bx, double* by, double* bz);

bool compareResults(int numTests, bool verbose = false);

int main() {
    const int numTests = 1000;
    bool result = compareResults(numTests, false);
    return result ? 0 : 1;
}

bool compareResults(int numTests, bool verbose) {
    std::cout << "Testing T89 model with " << numTests << " random samples..." << std::endl;

    bool pass = true;
    double tol = 2e-5;

    double maxDx = 0.0, maxDy = 0.0, maxDz = 0.0;

    std::srand(static_cast<unsigned>(std::time(nullptr)));

    for (int i = 0; i < numTests; ++i) {
        double x = -20.0 + static_cast<double>(rand()) / RAND_MAX * 30.0;
        double y = -10.0 + static_cast<double>(rand()) / RAND_MAX * 20.0;
        double z = -10.0 + static_cast<double>(rand()) / RAND_MAX * 20.0;
        double psiDeg = -35.0 + static_cast<double>(rand()) / RAND_MAX * 70.0;
        double psi = psiDeg * M_PI / 180.0;
        int iopt = 1 + rand() % 7;
        double parmod[10] = {0};  // unused

        double bx_f, by_f, bz_f;
        double bx_c, by_c, bz_c;

        t89c_(&iopt, parmod, &psi, &x, &y, &z, &bx_f, &by_f, &bz_f);
        t89(iopt, parmod, psi, x, y, z, &bx_c, &by_c, &bz_c);

        double dx = std::fabs(bx_f - bx_c);
        double dy = std::fabs(by_f - by_c);
        double dz = std::fabs(bz_f - bz_c);

        maxDx = std::max(maxDx, dx);
        maxDy = std::max(maxDy, dy);
        maxDz = std::max(maxDz, dz);

        if (dx > tol || dy > tol || dz > tol) {
            pass = false;
            if (verbose) {
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "Mismatch at test " << (i + 1) << "\n";
                std::cout << "  Input: x=" << x << " y=" << y << " z=" << z
                          << " psi=" << psi << " iopt=" << iopt << "\n";
                std::cout << "  Fortran:  Bx=" << bx_f << " By=" << by_f << " Bz=" << bz_f << "\n";
                std::cout << "  C++     : Bx=" << bx_c << " By=" << by_c << " Bz=" << bz_c << "\n";
                std::cout << "  Diff    : dx=" << dx << " dy=" << dy << " dz=" << dz << "\n\n";
            }
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Max differences: dx=" << maxDx
              << " dy=" << maxDy
              << " dz=" << maxDz << std::endl;

    if (pass) {
        std::cout << "All tests passed ✅" << std::endl;
    } else {
        std::cout << "Test failed ❌" << std::endl;
    }

    return pass;
}
