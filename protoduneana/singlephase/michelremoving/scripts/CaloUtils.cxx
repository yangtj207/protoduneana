#include "CaloUtils.h"
#include "TMath.h"

double beta(double gamma){
  double value=TMath::Sqrt(1-(1.0/(gamma*gamma)));
  return value;
}

double gamma(double KE,double mass){
  double value=(double(KE)/mass)+1;
  return value;
}

const double C=-5.2146;
const double X0=0.2;
const double X1=3.0;
const double a=0.19559;
const double m=3.0;
const double N=2*TMath::Log(10);

double density(double bg){//replaced x by x1
  double value;
  double x = TMath::Log10(bg);
  if(x<X0) return 0;
  if(x>X1) return N*x + C;
  value=a*(TMath::Power((X1-x),m));
  return N*x + C + value;
}

double Wmax(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double num=2*Me*(TMath::Power(b*g,2));
  double den=1+double(2*g*Me)/mass + TMath::Power((double(Me)/mass),2);
  double value=double(num)/den;
  return value;
}

double dpdx(double KE,double x,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double epsilon=(double(K)/2)*(double(Z)/A)*(double(x*rho)/(b*b));
  double A0=double(2*Me*(TMath::Power((b*g),2)))/I;
  double A1=double(epsilon)/I;
  double value=(1.0/x)*epsilon*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - b*b - density(b*g));
  return value;
}
