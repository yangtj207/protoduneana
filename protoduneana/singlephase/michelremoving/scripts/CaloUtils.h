#ifndef CALOUTILS
#define CALOUTILS

const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.511; // Mev for electron
const double rho=1.396;//g/cm^3

double beta(double gamma);

double gamma(double KE,double mass);

double density(double bg);

double Wmax(double KE,double mass);

double dpdx(double KE,double x,double mass);

double GetMuKEfromRange(double range);

double GetdEdx(double dQdx, double E_field, double calconst, double alp, double bet);

#endif
