#ifndef T96_MODEL_H
#define T96_MODEL_H

#include <array>
#include <vector>
#include <cmath>

// ---- Constants ----
constexpr double DTET0   = 0.034906;
constexpr double XLTDAY  = 78.0;
constexpr double XLTNGHT = 70.0;
constexpr double DEG2RAD = 0.01745329;
constexpr double PI      = 3.141592654;

// ---- Struct Declarations ----
struct Coord11;
struct RHDR;
struct LoopDip1;
struct Coord21;
struct DX1;

// ---- External Variables ----
extern Coord11 coord11;
extern LoopDip1 loopDip1;
extern RHDR rhdr;
extern Coord21 coord21;
extern DX1 dx1;

// ---- Core Function ----
void t96(int iopt, const std::array<double, 10>& parmod, double ps,
         double x, double y, double z, double& bx, double& by, double& bz);

// ---- Dipole Shielding ----
void t96DipoleShield(double psi, double x, double y, double z,
                     double *Bx, double *By, double *Bz);

void CylHarmPerp(double x, double y, double z,
                 double *a, double *b,
                 double *Bx, double *By, double *Bz);

void CylHarmPara(double x, double y, double z,
                 double *c, double *d,
                 double *Bx, double *By, double *Bz);

// ---- Interconnection Field ----
void t96Intercon(double x, double y, double z, double* Bx, double* By, double* Bz);

// ---- Tail Current Systems ----
void t96TailRC96(double sps, double x, double y, double z,
                 double& BXRC, double& BYRC, double& BZRC,
                 double& BXT2, double& BYT2, double& BZT2,
                 double& BXT3, double& BYT3, double& BZT3);

void t96RingCurr96(double x, double y, double z, 
                   double& Bx, double& By, double& Bz);

void t96TailDisk(double x, double y, double z,
                 double& Bx, double& By, double& Bz);

void t96Tail87(double x, double z, double* Bx, double* Bz);

void t96ShlCar3x3(const std::array<double, 48>& A, double x, double y, double z, 
                  double sps, double& HX, double& HY, double& HZ);

// ---- Region 1 Birkeland Currents ----
void t96Birk1Tot02(double ps, double x, double y, double z, double* Bx, double* By, double* Bz);

void t96DipLoop1(const std::array<double, 4>& xi, double D[3][26]);

void t96Circle(double x, double y, double z, double rl, 
               double* Bx, double* By, double* Bz);

void t96CrossLoop(double x, double y, double z, 
                  double* Bx, double* By, double* Bz, 
                  double xc, double rl, double al);

void t96DipXYZ(double x, double y, double z,
               double* Bxx, double* Byx, double* Bzx,
               double* Bxy, double* Byy, double* Bzy,
               double* Bxz, double* Byz, double* Bzz);

void t96Birk1Shield(double ps, double x, double y, double z, 
                    double* Bx, double* By, double* Bz);

// ---- Region 2 Birkeland Currents ----
void t96Birk2Tot_02(double ps, double x, double y, double z,
                    double* Bx, double* By, double* Bz);

void t96Birk2Shield(double x, double y, double z, double ps,
                    double* Hx, double* Hy, double* Hz);

void t96R2Birk(double x, double y, double z, double ps,
               double* Bx, double* By, double* Bz);

void t96R2Inner(double x, double y, double z, double* Bx, double* By, double* Bz);

void t96R2Outer(double x, double y, double z, double* Bx, double* By, double* Bz);

void t96R2Sheet(double x, double y, double z, double* Bx, double* By, double* Bz);

void t96Loops4(double x, double y, double z,
               double* Bx, double* By, double* Bz,
               double xc, double yc, double zc,
               double r, double theta, double phi);

// ---- Conical Harmonics & Dipole Distributions ----
void t96BConic(double x, double y, double z, 
               std::vector<double>& CBX, 
               std::vector<double>& CBY, 
               std::vector<double>& CBZ, 
               int NMAX);

void t96DipDistr(double x, double y, double z, double* Bx, double* By, double* Bz, int mode);

// ---- Utility Functions ----
double fexp(double ct, double p);   // Assuming you have this utility function
double xksi(double x, double y, double z);  // Assuming utility function for sheet transition
double tksi(double xks, double delarg, double delarg1);  // Transition function

#endif  // T96_MODEL_H
