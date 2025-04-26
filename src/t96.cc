#include "t96.h"

void t96(    int iopt, double *parmod, double psi, 
            double x, double y, double z,
            double *Bx, double *By, double *Bz) {

}


void t96DipoleShield(double psi, double x, double y, double z,
                        double *Bx, double *By, double *Bz) {
    
    /* calculate the magnetopause field which will
    shield the Earth's dipole field using cylindrical 
    harmonics. This bit comes from https://doi.org/10.1029/94JA03193 */

    /* the parameters to use for the cylindrical harmonics come
    from table 1 of the above paper*/
    double a[] = {0.24777,-27.003,-0.46815,7.0637,-1.5918,-0.090317};
    double b[] = {57.522,13.757,2.0100,10.458,4.5798,2.1695};
    double c[] = {-0.65385,-18.061,-0.40457,-5.0995,1.2846,.078231};
    double d[] = {39.592,13.291,1.9970,10.062,4.5140,2.1558};

    /* get the cylindrical harmonics for both 
    parallel and perpendicular components */
    double Bx0, By0, Bz0, Bx1, By1, Bz1;
    CylHarmPerp(x,y,z,a,b,&Bx0,&By0,&Bz0);
    CylHarmPara(x,y,z,a,b,&Bx1,&By1,&Bz1);


    /*combine them (equation 16)*/
    double cps = cos(psi);
    double sps = sin(psi);
    *Bx = Bx0*cps + Bx1*sps;
    *By = By0*cps + By1*sps;
    *Bz = Bz0*cps + Bz1*sps;

}

void CylHarmPerp(    double x, double y, double z,
                    double *a, double *b,
                    double *Bx, double *By, double *Bz) {
    
    /* I took this bit from the original Fortran code */
    double rho = sqrt(y*y + z*z);
    double sinp, cosp;
    if (rho < 1e-8) {
        sinp = 1.0;
        cosp = 0.0;
        rho = 1e-8;
    } else {
        sinp = z/rho;
        cosp = y/rho;
    }

    /* some variables which will be used more than once */
    double sinp2 = sinp*sinp;
    double cosp2 = cosp*cosp;
    double xb1,expxb0,expxb1,rhob0,rhob1,J0rb0,J0rb1,J1rb0,J1rb1;


    /* equation 10, 11 and 12 */
    double bx = 0.0, br = 0.0, bp = 0.0;
    bx = 0.0;
    br = 0.0;
    bp = 0.0;
    int i;
    for(i=0;i<3;i++) {
        /* get the common terms */
        xb1 = x/b[i+1];
        expxb0 = exp(x/b[i]);
        expxb1 = exp(xb1);
        rhob0 = rho/b[i];
        rhob1 = rho/b[i+3];
        J0rb0 = std::cyl_bessel_j(0,rhob0);
        J0rb1 = std::cyl_bessel_j(0,rhob1);
        J1rb0 = std::cyl_bessel_j(1,rhob0);
        J1rb1 = std::cyl_bessel_j(1,rhob1);

        /* sum them */
        bx += -a[i]*expxb0*J1rb0 + (a[i+3]/b[i+3])*expxb1*(rho*J0rb1 + x*J1rb1);
        br += a[i]*expxb0*(J1rb0/rhob0 - J1rb0) + a[i+3]*expxb1*(xb1*J0rb1 - (rhob1*rhob1 + xb1 - 1)*J1rb1/rhob1);
        bp += -a[i]*expxb0*J1rb0/rhob0 + a[i+3]*expxb1*(J0rb1 + ((x - b[i+3])/b[i+3])*J1rb1/rhob1);

    }
    /* multiply by sine or cosine */
    bx *= sinp;
    br *= sinp;
    bp *= cosp;

    /* convert back to GSM*/
    *Bx = bx;
    *By = br*cosp - bp*sinp;
    *Bz = br*sinp + bp*cosp;

}


void CylHarmPara(    double x, double y, double z,
                    double *c, double *d,
                    double *Bx, double *By, double *Bz) {
    
    /* I took this bit from the original Fortran code */
    double rho = sqrt(y*y + z*z);
    double sinp, cosp;
    if (rho < 1e-8) {
        sinp = 1.0;
        cosp = 0.0;
        rho = 1e-8;
    } else {
        sinp = z/rho;
        cosp = y/rho;
    }

    /* some variables which will be used more than once */
    double sinp2 = sinp*sinp;
    double cosp2 = cosp*cosp;
    double xd1,expxd0,expxd1,rhod0,rhod1,J0rd0,J0rd1,J1rd0,J1rd1;


    /* equation 13 and 14 (15 = 0)*/
    double bx = 0.0, br = 0.0, bp = 0.0;
    bx = 0.0;
    br = 0.0;

    int i;
    for(i=0;i<3;i++) {
        /* get the common terms */
        xd1 = x/d[i+1];
        expxd0 = exp(x/d[i]);
        expxd1 = exp(xd1);
        rhod0 = rho/d[i];
        rhod1 = rho/d[i+3];
        J0rd0 = std::cyl_bessel_j(0,rhod0);
        J0rd1 = std::cyl_bessel_j(0,rhod1);
        J1rd0 = std::cyl_bessel_j(1,rhod0);
        J1rd1 = std::cyl_bessel_j(1,rhod1);

        /* sum them */
        bx += -c[i]*expxd0*J1rd0 + c[i+3]*expxd1*(rhod1*J1rd1 -((x+d[i+3])/d[i+3])*J0rd1);
        br += -c[i]*expxd0*J1rd0 + (c[i+3]/d[i+3])*expxd1*(rho*J0rd1 - x*J1rd1);

    }

    /* convert back to GSM*/
    *Bx = bx;
    *By = br*cosp;
    *Bz = br*sinp;

}

void t96Intercon(    double x, double y, double z,
                    double *Bx, double *By, double *Bz) {


}

void t96RingCurrent() {

}

void t96TailDisk() {

}

void t96Tail87() {

}

void t96CartHarmonicShield() {

}

void t96Region1() {

}

void t96DipLoop() {

}

void t96Circle() {

}

void t96CrossLoop() {

}

void t96Dipolesxyz() {

}

void t96ConDip() {

}


/**
 * Computes the shielding magnetic field components (Bx, By, Bz)
 * for Region 1 Birkeland currents using box harmonics.
 *
 * @param ps   Dipole tilt angle in radians
 * @param x    GSM X coordinate (Re)
 * @param y    GSM Y coordinate (Re)
 * @param z    GSM Z coordinate (Re)
 * @param Bx   Pointer to store X component of magnetic field (nT)
 * @param By   Pointer to store Y component of magnetic field (nT)
 * @param Bz   Pointer to store Z component of magnetic field (nT)
 */
void t96Birk1Shield(double ps, double x, double y, double z, 
                    double* Bx, double* By, double* Bz) {
    // Initialize A array (80 coefficients)
    const std::array<double, 80> A = {{
        1.174198045, -1.463820502, 4.840161537, -3.674506864,
        82.18368896, -94.94071588, -4122.331796, 4670.278676,
        -21.54975037, 26.72661293, -72.81365728, 44.09887902,
        40.08073706, -51.23563510, 1955.348537, -1940.971550,
        794.0496433, -982.2441344, 1889.837171, -558.9779727,
        -1260.543238, 1260.063802, -293.5942373, 344.7250789,
        -773.7002492, 957.0094135, -1824.143669, 520.7994379,
        1192.484774, -1192.184565, 89.15537624, -98.52042999,
        -0.08168777675, 0.04255969908, 0.3155237661, -0.3841755213,
        2.494553332, -0.06571440817, -2.765661310, 0.4331001908,
        0.1099181537, -0.06154126980, -0.3258649260, 0.6698439193,
        -5.542735524, 0.1604203535, 5.854456934, -0.8323632049,
        3.732608869, -3.130002153, 107.0972607, -32.28483411,
        -115.2389298, 54.45064360, -0.5826853320, -3.582482231,
        -4.046544561, 3.311978102, -104.0839563, 30.26401293,
        97.29109008, -50.62370872, -296.3734955, 127.7872523,
        5.303648988, 10.40368955, 69.65230348, 466.5099509,
        1.645049286, 3.825838190, 11.66675599, 558.9781177,
        1.826531343, 2.066018073, 25.40971369, 990.2795225,
        2.319489258, 4.555148484, 9.691185703, 591.8280358
    }};

    // P1, R1, Q1, S1 mapped via equivalence to A[64] to A[79]
    const double* P1 = &A[64];
    const double* R1 = &A[68];
    const double* Q1 = &A[72];
    const double* S1 = &A[76];

    double RP[4], RR[4], RQ[4], RS[4];

    double cps = std::cos(ps);
    double sps = std::sin(ps);
    double s3ps = 4.0 * cps * cps - 1.0;

    *Bx = 0.0;
    *By = 0.0;
    *Bz = 0.0;

    // Precompute reciprocals
    for (int i = 0; i < 4; ++i) {
        RP[i] = 1.0 / P1[i];
        RR[i] = 1.0 / R1[i];
        RQ[i] = 1.0 / Q1[i];
        RS[i] = 1.0 / S1[i];
    }

    int L = 0;

    for (int m = 1; m <= 2; ++m) {  // m = symmetry type
        for (int i = 0; i < 4; ++i) {
            double cypi = std::cos(y * RP[i]);
            double cyqi = std::cos(y * RQ[i]);
            double sypi = std::sin(y * RP[i]);
            double syqi = std::sin(y * RQ[i]);

            for (int k = 0; k < 4; ++k) {
                double szrk = std::sin(z * RR[k]);
                double czsk = std::cos(z * RS[k]);
                double czrk = std::cos(z * RR[k]);
                double szsk = std::sin(z * RS[k]);

                double sqpr = std::sqrt(RP[i]*RP[i] + RR[k]*RR[k]);
                double sqqs = std::sqrt(RQ[i]*RQ[i] + RS[k]*RS[k]);

                double epr = std::exp(x * sqpr);
                double eqs = std::exp(x * sqqs);

                for (int n = 1; n <= 2; ++n) {
                    double hx, hy, hz;

                    if (m == 1) {  // Perpendicular symmetry
                        if (n == 1) {
                            hx = -sqpr * epr * cypi * szrk;
                            hy = RP[i] * epr * sypi * szrk;
                            hz = -RR[k] * epr * cypi * czrk;
                        } else {
                            hx *= cps;
                            hy *= cps;
                            hz *= cps;
                        }
                    } else {  // Parallel symmetry
                        if (n == 1) {
                            hx = -sps * sqqs * eqs * cyqi * czsk;
                            hy = sps * RQ[i] * eqs * syqi * czsk;
                            hz = sps * RS[k] * eqs * cyqi * szsk;
                        } else {
                            hx *= s3ps;
                            hy *= s3ps;
                            hz *= s3ps;
                        }
                    }

                    *Bx += A[L] * hx;
                    *By += A[L] * hy;
                    *Bz += A[L] * hz;
                    ++L;
                }
            }
        }
    }
}

void t96Region2() {

}

void t96Region2Shield() {

}

/**
 * Computes the magnetic field components for the inner part of Region 2 field-aligned currents.
 */
void t96R2Inner(double x, double y, double z, double* Bx, double* By, double* Bz) {
    // Coefficients from Fortran DATA statements
    const double PL[8] = {154.185, -2.12446, 0.0601735, -0.00153954, 0.0000355077, 
                          29.9996, 262.886, 99.9132};

    const double PN[8] = {-8.1902, 6.5239, 5.504, 7.7815, 0.8573, 3.0986, 0.0774, -0.038};

    // Arrays for conical harmonics
    std::vector<double> cbx(5), cby(5), cbz(5);

    // Compute conical harmonics
    t96BConic(x, y, z, cbx, cby, cbz, 5);

    // Compute 4-loop current system contribution
    double dbx8, dby8, dbz8;
    t96Loops4(x, y, z, &dbx8, &dby8, &dbz8, PN[0], PN[1], PN[2], PN[3], PN[4], PN[5]);

    // Compute dipole distribution contributions
    double dbx6, dby6, dbz6;
    double dbx7, dby7, dbz7;

    t96DipDistr(x - PN[6], y, z, &dbx6, &dby6, &dbz6, 0);
    t96DipDistr(x - PN[7], y, z, &dbx7, &dby7, &dbz7, 1);

    // Combine all contributions to get total field components
    *Bx = PL[0]*cbx[0] + PL[1]*cbx[1] + PL[2]*cbx[2] + PL[3]*cbx[3] + PL[4]*cbx[4]
        + PL[5]*dbx6 + PL[6]*dbx7 + PL[7]*dbx8;

    *By = PL[0]*cby[0] + PL[1]*cby[1] + PL[2]*cby[2] + PL[3]*cby[3] + PL[4]*cby[4]
        + PL[5]*dby6 + PL[6]*dby7 + PL[7]*dby8;

    *Bz = PL[0]*cbz[0] + PL[1]*cbz[1] + PL[2]*cbz[2] + PL[3]*cbz[3] + PL[4]*cbz[4]
        + PL[5]*dbz6 + PL[6]*dbz7 + PL[7]*dbz8;
}


/**
 * Computes conical harmonics magnetic field components (CBX, CBY, CBZ)
 * up to order NMAX at position (x, y, z) in GSM coordinates.
 *
 * @param x     X coordinate (Re)
 * @param y     Y coordinate (Re)
 * @param z     Z coordinate (Re)
 * @param CBX   Vector to store X components [size NMAX]
 * @param CBY   Vector to store Y components [size NMAX]
 * @param CBZ   Vector to store Z components [size NMAX]
 * @param NMAX  Number of harmonic terms
 */
void t96BConic(double x, double y, double z, 
               std::vector<double>& CBX, 
               std::vector<double>& CBY, 
               std::vector<double>& CBZ, 
               int NMAX) {
    
    double ro2 = x * x + y * y;
    double ro = std::sqrt(ro2);

    if (ro == 0.0) {
        // Avoid division by zero at the Z-axis
        CBX.assign(NMAX, 0.0);
        CBY.assign(NMAX, 0.0);
        CBZ.assign(NMAX, 0.0);
        return;
    }

    double cf = x / ro;
    double sf = y / ro;

    double cfm1 = 1.0;
    double sfm1 = 0.0;

    double r2 = ro2 + z * z;
    double r = std::sqrt(r2);

    double c = z / r;
    double s = ro / r;

    double ch = std::sqrt(0.5 * (1.0 + c));
    double sh = std::sqrt(0.5 * (1.0 - c));

    double tnhm1 = 1.0;
    double cnhm1 = 1.0;

    double tnh = sh / ch;
    double cnh = 1.0 / tnh;

    // Resize vectors if needed
    if (CBX.size() != NMAX) CBX.resize(NMAX);
    if (CBY.size() != NMAX) CBY.resize(NMAX);
    if (CBZ.size() != NMAX) CBZ.resize(NMAX);

    for (int m = 1; m <= NMAX; ++m) {
        double cfm = cfm1 * cf - sfm1 * sf;
        double sfm = cfm1 * sf + sfm1 * cf;

        cfm1 = cfm;
        sfm1 = sfm;

        double tnhm = tnhm1 * tnh;
        double cnhm = cnhm1 * cnh;

        double bt = m * cfm / (r * s) * (tnhm + cnhm);
        double bf = -0.5 * m * sfm / r * (tnhm1 / (ch * ch) - cnhm1 / (sh * sh));

        tnhm1 = tnhm;
        cnhm1 = cnhm;

        CBX[m - 1] = bt * c * cf - bf * sf;
        CBY[m - 1] = bt * c * sf + bf * cf;
        CBZ[m - 1] = -bt * s;
    }
}


/**
 * Computes the magnetic field components (Bx, By, Bz) from a linear distribution 
 * of dipolar sources along the Z-axis, as used in the Tsyganenko 1996 (T96) model.
 *
 * @param x   X coordinate in GSM (Re)
 * @param y   Y coordinate in GSM (Re)
 * @param z   Z coordinate in GSM (Re)
 * @param mode  Determines dipole strength variation:
 *              0 = Step-function distribution
 *              1 = Linear variation
 * @param Bx  Pointer to store X component of the magnetic field (nT)
 * @param By  Pointer to store Y component of the magnetic field (nT)
 * @param Bz  Pointer to store Z component of the magnetic field (nT)
 */
void t96DipDistr(double x, double y, double z, double* Bx, double* By, double* Bz, int mode) {
    double x2 = x * x;
    double rho2 = x2 + y * y;
    double r2 = rho2 + z * z;

    if (rho2 == 0.0) {
        // Avoid division by zero if point is exactly on Z-axis
        *Bx = 0.0;
        *By = 0.0;
        *Bz = 0.0;
        return;
    }

    if (mode == 0) {
        double r3 = r2 * std::sqrt(r2);
        *Bx = z / (rho2 * rho2) * (r2 * (y * y - x2) - rho2 * x2) / r3;
        *By = -x * y * z / (rho2 * rho2) * (2.0 * r2 + rho2) / r3;
        *Bz = x / r3;
    } else {
        *Bx = z / (rho2 * rho2) * (y * y - x2);
        *By = -2.0 * x * y * z / (rho2 * rho2);
        *Bz = x / rho2;
    }
}

void t96Region2Outer() {

}

void t964CurrentLoops() {

}

void t96R2Sheet(double x, double y, double z, double* Bx, double* By, double* Bz) {
    static const double PNONX[] = {
        19.0969, -9.28828, -0.129687, 5.58594,
        22.5055, 0.48375e-01, 0.396953e-01, 0.579023e-01
    };
    static const double PNONY[] = {
        -13.6750, -6.70625, 2.31875, 11.4062,
        20.4562, 0.478750e-01, 0.363750e-01, 0.567500e-01
    };
    static const double PNONZ[] = {
        -16.7125, -16.4625, -0.1625, 5.1, 
        23.7125, 0.355625e-01, 0.318750e-01, 0.53875e-01
    };
    
    const double A[80] = {
        8.07190, -7.39582, -7.62341, 0.684671, -13.5672, 11.6681, 13.1154, -0.890217,
        7.78726, -5.38346, -8.08738, 0.609385, -2.70410, 3.53741, 3.15549, -1.11069,
       -8.47555, 0.278122, 2.73514, 4.55625, 13.1134, 1.15848, -3.52648, -8.24698,
       -6.85710, -2.81369, 2.03795, 4.64383, 2.49309, -1.22041, -1.67432, -0.422526,
       -5.39796, 7.10326, 5.53730, -13.1918, 4.67853, -7.60329, -2.53066, 7.76338,
        5.60165, 5.34816, -4.56441, 7.05976, -2.62723, -0.529078, 1.42019, -2.93919,
       55.6338, -1.55181, 39.8311, -80.6561, -46.9655, 32.8925, -6.32296, 19.7841,
       124.731, 10.4347, -30.7581, 102.680, -47.4037, -3.31278, 9.37141, -50.0268,
      -533.319, 110.426, 1000.20, -1051.40, 1619.48, 589.855, -1462.73, 1087.10,
      -1994.73, -1654.12, 1263.33, -260.210, 1424.84, 1255.71, -956.733, 219.946
   };
   
   const double B[80] = {
       -9.08427, 10.6777, 10.3288, -0.969987, 6.45257, -8.42508, -7.97464, 1.41996,
       -1.92490, 3.93575, 2.83283, -1.48621, 0.244033, -0.757941, -0.386557, 0.344566,
        9.56674, -2.5365, -3.32916, -5.86712, -6.19625, 1.83879, 2.52772, 4.34417,
        1.87268, -2.13213, -1.69134, -0.176379, -0.261359, 0.566419, 0.3138, -0.134699,
       -3.83086, -8.4154, 4.77005, -9.31479, 37.5715, 19.3992, -17.9582, 36.4604,
      -14.9993, -3.1442, 6.17409, -15.5519, 2.28621, -0.00891549, -0.462912, 2.47314,
       41.7555, 208.614, -45.7861, -77.8687, 239.357, -67.9226, 66.8743, 238.534,
      -112.136, 16.2069, -40.4706, -134.328, 21.56, -0.201725, 2.21, 32.5855,
      -108.217, -1005.98, 585.753, 323.668, -817.056, 235.750, -560.965, -576.892,
       684.193, 85.0275, 168.394, 477.776, -289.253, -123.216, 75.6501, -178.605
   };
   
   const double C[80] = {
      1167.61, -917.782, -1253.2, -274.128, -1538.75, 1257.62, 1745.07, 113.479,
       393.326, -426.858, -641.1, 190.833, -29.9435, -1.04881, 117.125, -25.7663,
      -1168.16, 910.247, 1239.31, 289.515, 1540.56, -1248.29, -1727.61, -131.785,
      -394.577, 426.163, 637.422, -187.965, 30.0348, 0.221898, -116.68, 26.0291,
        12.6804, 4.84091, 1.18166, -2.75946, -17.9822, -6.80357, -1.47134, 3.02266,
         4.79648, 0.665255, -0.256229, -0.0857282, -0.588997, 0.0634812, 0.164303, -0.15285,
        22.2524, -22.4376, -3.85595, 6.07625, -105.959, -41.6698, 0.378615, 1.55958,
        44.3981, 18.8521, 3.19466, 5.89142, -8.63227, -2.36418, -1.027, -2.31515,
      1035.38, 2040.66, -131.881, -744.533, -3274.93, -4845.61, 482.438, 1567.43,
      1354.02, 2040.47, -151.653, -845.012, -111.723, -265.343, -26.1171, 216.632
   };


   
    double xks = xksi(x, y, z);  // Variation across the current sheet

    // Compute T1, T2, T3 terms for X, Y, Z
    double t1x = xks / std::sqrt(xks * xks + PNONX[5] * PNONX[5]);
    double t2x = std::pow(PNONX[6], 3) / std::pow(std::sqrt(xks * xks + PNONX[6] * PNONX[6]), 3);
    double t3x = xks / std::pow(std::sqrt(xks * xks + PNONX[7] * PNONX[7]), 5) * 3.493856 * std::pow(PNONX[7], 4);

    double t1y = xks / std::sqrt(xks * xks + PNONY[5] * PNONY[5]);
    double t2y = std::pow(PNONY[6], 3) / std::pow(std::sqrt(xks * xks + PNONY[6] * PNONY[6]), 3);
    double t3y = xks / std::pow(std::sqrt(xks * xks + PNONY[7] * PNONY[7]), 5) * 3.493856 * std::pow(PNONY[7], 4);

    double t1z = xks / std::sqrt(xks * xks + PNONZ[5] * PNONZ[5]);
    double t2z = std::pow(PNONZ[6], 3) / std::pow(std::sqrt(xks * xks + PNONZ[6] * PNONZ[6]), 3);
    double t3z = xks / std::pow(std::sqrt(xks * xks + PNONZ[7] * PNONZ[7]), 5) * 3.493856 * std::pow(PNONZ[7], 4);

    // Cylindrical and spherical coordinates
    double rho2 = x * x + y * y;
    double r = std::sqrt(rho2 + z * z);
    double rho = std::sqrt(rho2);

    // Azimuthal angle terms
    double c1p = x / rho;
    double s1p = y / rho;
    double s2p = 2.0 * s1p * c1p;
    double c2p = c1p * c1p - s1p * s1p;
    double s3p = s2p * c1p + c2p * s1p;
    double c3p = c2p * c1p - s2p * s1p;
    double s4p = s3p * c1p + c3p * s1p;

    double ct = z / r;
    double st = rho / r;

    // Compute S1 to S5 for BX using PNONX
    double s1 = fexp(ct, PNONX[0]);
    double s2 = fexp(ct, PNONX[1]);
    double s3 = fexp(ct, PNONX[2]);
    double s4 = fexp(ct, PNONX[3]);
    double s5 = fexp(ct, PNONX[4]);

    // Compute BX component
    *Bx = 
        s1 * ((A[0] + A[1]*t1x + A[2]*t2x + A[3]*t3x)
            + c1p * (A[4] + A[5]*t1x + A[6]*t2x + A[7]*t3x)
            + c2p * (A[8] + A[9]*t1x + A[10]*t2x + A[11]*t3x)
            + c3p * (A[12] + A[13]*t1x + A[14]*t2x + A[15]*t3x))
        + s2 * ((A[16] + A[17]*t1x + A[18]*t2x + A[19]*t3x)
            + c1p * (A[20] + A[21]*t1x + A[22]*t2x + A[23]*t3x)
            + c2p * (A[24] + A[25]*t1x + A[26]*t2x + A[27]*t3x)
            + c3p * (A[28] + A[29]*t1x + A[30]*t2x + A[31]*t3x))
        + s3 * ((A[32] + A[33]*t1x + A[34]*t2x + A[35]*t3x)
            + c1p * (A[36] + A[37]*t1x + A[38]*t2x + A[39]*t3x)
            + c2p * (A[40] + A[41]*t1x + A[42]*t2x + A[43]*t3x)
            + c3p * (A[44] + A[45]*t1x + A[46]*t2x + A[47]*t3x))
        + s4 * ((A[48] + A[49]*t1x + A[50]*t2x + A[51]*t3x)
            + c1p * (A[52] + A[53]*t1x + A[54]*t2x + A[55]*t3x)
            + c2p * (A[56] + A[57]*t1x + A[58]*t2x + A[59]*t3x)
            + c3p * (A[60] + A[61]*t1x + A[62]*t2x + A[63]*t3x))
        + s5 * ((A[64] + A[65]*t1x + A[66]*t2x + A[67]*t3x)
            + c1p * (A[68] + A[69]*t1x + A[70]*t2x + A[71]*t3x)
            + c2p * (A[72] + A[73]*t1x + A[74]*t2x + A[75]*t3x)
            + c3p * (A[76] + A[77]*t1x + A[78]*t2x + A[79]*t3x));

    // Repeat for BY using PNONY and B array
    s1 = fexp(ct, PNONY[0]);
    s2 = fexp(ct, PNONY[1]);
    s3 = fexp(ct, PNONY[2]);
    s4 = fexp(ct, PNONY[3]);
    s5 = fexp(ct, PNONY[4]);

    *By = 
        s1 * (s1p*(B[0]+B[1]*t1y+B[2]*t2y+B[3]*t3y)
            + s2p*(B[4]+B[5]*t1y+B[6]*t2y+B[7]*t3y)
            + s3p*(B[8]+B[9]*t1y+B[10]*t2y+B[11]*t3y)
            + s4p*(B[12]+B[13]*t1y+B[14]*t2y+B[15]*t3y))
        + s2 * (s1p*(B[16]+B[17]*t1y+B[18]*t2y+B[19]*t3y)
            + s2p*(B[20]+B[21]*t1y+B[22]*t2y+B[23]*t3y)
            + s3p*(B[24]+B[25]*t1y+B[26]*t2y+B[27]*t3y)
            + s4p*(B[28]+B[29]*t1y+B[30]*t2y+B[31]*t3y))
        + s3 * (s1p*(B[32]+B[33]*t1y+B[34]*t2y+B[35]*t3y)
            + s2p*(B[36]+B[37]*t1y+B[38]*t2y+B[39]*t3y)
            + s3p*(B[40]+B[41]*t1y+B[42]*t2y+B[43]*t3y)
            + s4p*(B[44]+B[45]*t1y+B[46]*t2y+B[47]*t3y))
        + s4 * (s1p*(B[48]+B[49]*t1y+B[50]*t2y+B[51]*t3y)
            + s2p*(B[52]+B[53]*t1y+B[54]*t2y+B[55]*t3y)
            + s3p*(B[56]+B[57]*t1y+B[58]*t2y+B[59]*t3y)
            + s4p*(B[60]+B[61]*t1y+B[62]*t2y+B[63]*t3y))
        + s5 * (s1p*(B[64]+B[65]*t1y+B[66]*t2y+B[67]*t3y)
            + s2p*(B[68]+B[69]*t1y+B[70]*t2y+B[71]*t3y)
            + s3p*(B[72]+B[73]*t1y+B[74]*t2y+B[75]*t3y)
            + s4p*(B[76]+B[77]*t1y+B[78]*t2y+B[79]*t3y));

    // Finally, BZ using PNONZ and C array with fexp1
    s1 = fexp1(ct, PNONZ[0]);
    s2 = fexp1(ct, PNONZ[1]);
    s3 = fexp1(ct, PNONZ[2]);
    s4 = fexp1(ct, PNONZ[3]);
    s5 = fexp1(ct, PNONZ[4]);

    *Bz = 
        s1 * ((C[0]+C[1]*t1z+C[2]*t2z+C[3]*t3z)
            + c1p*(C[4]+C[5]*t1z+C[6]*t2z+C[7]*t3z)
            + c2p*(C[8]+C[9]*t1z+C[10]*t2z+C[11]*t3z)
            + c3p*(C[12]+C[13]*t1z+C[14]*t2z+C[15]*t3z))
        + s2 * ((C[16]+C[17]*t1z+C[18]*t2z+C[19]*t3z)
            + c1p*(C[20]+C[21]*t1z+C[22]*t2z+C[23]*t3z)
            + c2p*(C[24]+C[25]*t1z+C[26]*t2z+C[27]*t3z)
            + c3p*(C[28]+C[29]*t1z+C[30]*t2z+C[31]*t3z))
        + s3 * ((C[32]+C[33]*t1z+C[34]*t2z+C[35]*t3z)
            + c1p*(C[36]+C[37]*t1z+C[38]*t2z+C[39]*t3z)
            + c2p*(C[40]+C[41]*t1z+C[42]*t2z+C[43]*t3z)
            + c3p*(C[44]+C[45]*t1z+C[46]*t2z+C[47]*t3z))
        + s4 * ((C[48]+C[49]*t1z+C[50]*t2z+C[51]*t3z)
            + c1p*(C[52]+C[53]*t1z+C[54]*t2z+C[55]*t3z)
            + c2p*(C[56]+C[57]*t1z+C[58]*t2z+C[59]*t3z)
            + c3p*(C[60]+C[61]*t1z+C[62]*t2z+C[63]*t3z))
        + s5 * ((C[64]+C[65]*t1z+C[66]*t2z+C[67]*t3z)
            + c1p*(C[68]+C[69]*t1z+C[70]*t2z+C[71]*t3z)
            + c2p*(C[72]+C[73]*t1z+C[74]*t2z+C[75]*t3z)
            + c3p*(C[76]+C[77]*t1z+C[78]*t2z+C[79]*t3z));
}
   

/*
 * Function: xksi
 * ------------------------
 * Computes a scalar coordinate transformation used in the Tsyganenko 1996 (T96)
 * geomagnetic field model.
 *
 * This function is part of the "stretching" system of coordinates that maps
 * dipole field lines into a stretched tail-like geometry, better reflecting
 * the real shape of the magnetosphere.
 *
 * Inputs:
 *   x, y, z : GSM position coordinates in Earth radii (Re)
 *
 * Returns:
 *   XKSI : a transformed scalar that incorporates stretching of coordinates 
 *          based on the position vector and empirical parameters. It reflects
 *          the difference between a geometrical scaling factor (α) and an angular
 *          factor (ϕ) related to latitude within the stretched geometry.
 *
 * References:
 *   Tsyganenko, N. A. (1995, 1996) - NB#3 (notebooks), P.26-27
 */
double xksi(double x, double y, double z) {
    // Empirical stretch coefficients (Axx, Bxx, Cxx, R0, DR)
    const double A11A12 = 0.305662, A21A22 = -0.383593, A41A42 = 0.2677733;
    const double A51A52 = -0.097656, A61A62 = -0.636034;
    const double B11B12 = -0.359862, B21B22 = 0.424706;
    const double C61C62 = -0.126366, C71C72 = 0.292578;
    const double R0 = 1.21563, DR = 7.50937;

    // Noon-midnight latitude mapping (in radians)
    const double TNOON = 0.3665191;   // ~69 degrees
    const double DTETA = 0.09599309;  // width ~5.5 degrees

    // Precompute DR² to avoid recomputation
    double dr2 = DR * DR;

    // Compute radial distance and its powers
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double r2 = x2 + y2 + z2;
    double r = std::sqrt(r2);

    // Normalize position vector
    double xr = x / r;
    double yr = y / r;
    double zr = z / r;

    // Compute position-dependent stretch magnitude (PR)
    // PR increases with distance beyond R0 in a smoothed way
    double pr;
    if (r < R0) {
        pr = 0.0;
    } else {
        pr = std::sqrt((r - R0) * (r - R0) + dr2) - DR;
    }

    // Stretched coordinates (F, G, H) — warped versions of x, y, z
    double f = x + pr * (A11A12 + A21A22 * xr + A41A42 * xr * xr +
                         A51A52 * yr * yr + A61A62 * zr * zr);
    double g = y + pr * (B11B12 * yr + B21B22 * xr * yr);
    double h = z + pr * (C61C62 * zr + C71C72 * xr * zr);

    // Prepare composite metrics from F, G, H
    double f2 = f * f;
    double g2 = g * g;
    double h2 = h * h;

    double fgh = f2 + g2 + h2;                   // squared length of (f, g, h)
    double fgh32 = std::pow(std::sqrt(fgh), 3);  // |FGH|^3

    double fchsq = f2 + g2;                      // projection onto the "equatorial" plane

    // Special case: avoid division by zero on the Z-axis (f² + g² ~ 0)
    if (fchsq < 1e-5) {
        return -1.0;
    }

    // Compute stretched parameters
    double sqfchsq = std::sqrt(fchsq);
    double alpha = fchsq / fgh32;  // geometric scaling factor

    // Map F coordinate into an effective magnetic latitude angle
    double theta = TNOON + 0.5 * DTETA * (1.0 - f / sqfchsq);
    double phi = std::pow(std::sin(theta), 2);  // angular mapping term

    // Final scalar transformation
    return alpha - phi;  // this value is often used as an argument in tapering functions
}

double fexp(double s, double a) {
    const double E = 2.718281828459;

    if (a < 0.0) {
        return std::sqrt(-2.0 * a * E) * s * std::exp(a * s * s);
    } else {
        return s * std::exp(a * (s * s - 1.0));
    }
}

double fexp1(double s, double a) {
    if (a <= 0.0) {
        return std::exp(a * s * s);
    } else {
        return std::exp(a * (s * s - 1.0));
    }
}

double tksi(double xksi, double xks0, double dxksi) {
    static bool initialized = false;
    static double tdz3;

    if (!initialized) {
        tdz3 = 2.0 * std::pow(dxksi, 3);
        initialized = true;
    }

    double tksii = 0.0;

    if (xksi - xks0 < -dxksi) {
        tksii = 0.0;
    } else if (xksi - xks0 >= dxksi) {
        tksii = 1.0;
    } else if (xksi >= xks0 - dxksi && xksi < xks0) {
        double br3 = std::pow(xksi - xks0 + dxksi, 3);
        tksii = 1.5 * br3 / (tdz3 + br3);
    } else if (xksi >= xks0 && xksi < xks0 + dxksi) {
        double br3 = std::pow(xksi - xks0 - dxksi, 3);
        tksii = 1.0 + 1.5 * br3 / (tdz3 - br3);
    }

    return tksii;
}



void t96Dipole(double psi, double x, double y, double z, 
               double* Bx, double* By, double* Bz) {
    
    // Precompute trigonometric values
    double sps = std::sin(psi);
    double cps = std::cos(psi);

    // Intermediate variables
    double P = x * x;
    double T = y * y;
    double U = z * z;
    double V = 3.0 * z * x;

    // Compute dipole scaling factor (constant from Fortran)
    double R2 = P + T + U;
    double Q = 30574.0 / std::pow(R2, 2.5);  // Equivalent to sqrt(R2)^5

    // Compute components in GSM
    *Bx = Q * ((T + U - 2.0 * P) * sps - V * cps);
    *By = -3.0 * y * Q * (x * sps + z * cps);
    *Bz = Q * ((P + T - 2.0 * U) * cps - V * sps);
}
