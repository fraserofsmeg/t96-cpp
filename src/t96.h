#include <cmath>
#include <vector>
#include <array>


struct WarpData {
    double CPSS, SPSS, DPSRR, RPS, WARP, D;
    double XS, ZS, DXSX, DXSY, DXSZ;
    double DZSX, DZSY, DZSZ;
    double DZETAS, DDZETADX, DDZETADY, DDZETADZ, ZSWW;
};

// Create a global instance (similar to Fortran COMMON block)
WarpData warpData;
