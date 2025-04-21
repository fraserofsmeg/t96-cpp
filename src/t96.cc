#include "t96.h"

void t96(	int iopt, double *parmod, double psi, 
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
	*Bx = Bx0*cps + Bx0*sps;
	*By = By0*cps + By0*sps;
	*Bz = Bz0*cps + Bz0*sps;

}

void CylHarmPerp(	double x, double y, double z,
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
		J0rb0 = j0(rhob0);
		J0rb1 = j0(rhob1);
		J1rb0 = j1(rhob0);
		J1rb1 = j1(rhob1);

		/* sum them */
		bx += -a[i]*expxb0*J1rb0 + (a[i+3]/b[i+3])*expxb1*(rho*J0rb1 + x*J1rb1);
		br += a[i]*expxb0*(J1rb0/rhob0 - J1rb0) + a[i+3]*expxb1*(xb1*J0rb1 - (rhob1*rhob1 + xb1 - 1)*J1rb1/rhob1);
		bp += -a[i]*expxb0*J1rb0/rhob0 + a[i+3]*expxb1*(J0rb1 + ((x - b[i+3])/b[i+3])*J1rb1/rhob1);

	}
	/* multiply by sine or cosine */
	bx *= sinp;
	br *= sinp;
	Bp *= cosp;

	/* convert back to GSM*/
	*Bx = bx;
	*By = br*cosp - bp*sinp;
	*Bz = br*sinp + bp*cosp;

}


void CylHarmPara(	double x, double y, double z,
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
		J0rd0 = j0(rhod0);
		J0rd1 = j0(rhod1);
		J1rd0 = j1(rhod0);
		J1rd1 = j1(rhod1);

		/* sum them */
		bx += -c[i]*expxd0*J1rd0 + c[i+3]*expxd1*(rhod1*J1rd1 -((x+d[i+3])/d[i+3])*J0rd1);
		br += -c[i]*expxd0*J1rd0 + (c[i+3]/d[i+3])*expxd1*(rho*J0rd1 - x*J1rd1);

	}

	/* convert back to GSM*/
	*Bx = bx;
	*By = br*cosp;
	*Bz = br*sinp;

}

void t96Intercon(	double x, double y, double z,
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

void t96Region1Shield() {

}

void t96Region2() {

}

void t96Region2Shield() {

}

void t96Region2Inner() {

}

void t96ConicalHarmonics() {

}

void t96DipoleDist() {

}

void t96Region2Outer() {

}

void t964CurrentLoops() {

}

void t96Region2Sheet() {

}

double xksi() {

}

double tksi() {

}

void t96Dipole( double psi, double x, double y, double z, 
				double *Bx, souble *By, double *Bz) {
	
	
}