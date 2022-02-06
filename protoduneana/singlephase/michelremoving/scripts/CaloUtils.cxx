#include "CaloUtils.h"
#include "TMath.h"
#include <TSpline.h>

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

double GetMuKEfromRange(double range){

  const int np = 13;
  double spline_KE[np] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000};
  double spline_Range[np] = {0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476};

  TSpline3 *sp = new TSpline3("Cubic Spline", spline_Range, spline_KE,13,"b2e2",0,0);

  double KE = sp->Eval(range);
  
  delete sp;
  return KE;
}

double GetdEdx(double dQdx, double E_field, double calconst, double alp, double bet){

  double Beta = bet/(rho*E_field);
  double Wion = 23.6e-6;
  return (exp(Beta * Wion *dQdx/calconst) - alp) / Beta;

}
