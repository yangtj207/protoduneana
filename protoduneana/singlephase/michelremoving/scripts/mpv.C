const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double rho=1.396;//g/cm^3
const double Mmu =105.6583755; //MeV
const double Me = 0.51099895000; //MeV
const double x = 0.5; //cm

const double C=-5.2146;
const double X0=0.2;
const double X1=3.0;
const double a=0.19559;
const double m=3.0;
const double N=2*TMath::Log(10);

double density(double bg){
  double value;
  double x = TMath::Log10(bg);
  if(x<X0) return 0;
  if(x>X1) return N*x + C;
  value=a*(TMath::Power((X1-x),m));
  return N*x + C + value;
}

ROOT::Math::VavilovAccurate vav;

void mpv(){

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  
  //https://pdg.lbl.gov/2011/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
  const int n1 = 15;
  
  double vp[15] = { 4.704e1,
                    5.616e1,
                    6.802e1,
                    8.509e1,
                    1.003e2,
                    1.527e2,
                    1.764e2,
                    2.218e2,
                    2.868e2,
                    3.567e2,
                    3.917e2,
                    4.945e2,
                    8.995e2,
                    1.101e3,
                    1.502e3};
  
  double vdedx[15] = {5.687,
                      4.461,
                      3.502,
                      2.731,
                      2.340,
                      1.771,
                      1.669,
                      1.570,
                      1.518,
                      1.508,
                      1.509,
                      1.526,
                      1.610,
                      1.644,
                      1.699};

  TGraph *grdedx = new TGraph(n1, vp, vdedx);
//  grdedx->Draw("ap");
//  gPad->SetLogx();
//  gPad->SetLogy();

  //const int n2 = 30;
  const int n2 = 100;
  double KE[n2];
  double pdgdEdx[n2];
  double VavdEdx[n2];
  double vkappa[n2];
  for (int i = 0; i<n2; ++i){
    KE[i] = 5*(i+1);
    double Em = KE[i] + Mmu;
    double Pm = sqrt(Em*Em - Mmu*Mmu);
    double gamma = Em/Mmu;
    double beta = sqrt(1-1/(gamma*gamma));
    double Xi=(double(K)/2)*(double(Z)/A)*(double(x*rho)/(beta*beta));
    double avgde = grdedx->Eval(Pm)*rho*x;
    double Emax = 2*Me*beta*beta*gamma*gamma/(1+2*gamma*Me/Mmu+(Me*Me)/(Mmu*Mmu));
    double kappa = Xi/Emax;
    vkappa[i] = kappa;
    //cout<<"KE = "<<KE[i]<<" MeV, Emax = "<<Emax<<" MeV, Xi = "<<Xi<<" MeV, kappa = "<<kappa<<", dE = "<<avgde<<" MeV"<<endl;
    //double Lambda = (de-avgde)/Xi - 0.422784 - beta*beta - log(Xi/Emax);
    //vav.Mode() returns the value of λ where the pdf is maximal function, and set kappa and beta2, if necessary.
    double mpvde = (vav.Mode(kappa, beta*beta) + 0.422784 + beta*beta + log(Xi/Emax))*Xi+avgde;
 
    //Pdg calculation
    double A0=double(2*Me*(TMath::Power((beta*gamma),2)))/I;
    double A1=double(Xi)/I;
    double pdgde = Xi*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - beta*beta - density(beta*gamma));
    //cout<<"KE = "<<KE[i]<<" MeV, Vavilov MPV dE/dx = "<<mpvde/x<<" MeV/cm, PDG MPV dEdx = "<<pdgde/x<<endl;
    pdgdEdx[i] = pdgde/x;
    VavdEdx[i] = mpvde/x;
    
  }

  TGraph *grPdgdEdx = new TGraph(n2, KE, pdgdEdx);
  TGraph *grVavdEdx = new TGraph(n2, KE, VavdEdx);

  TCanvas *c1 = new TCanvas("c1","c1");
  grPdgdEdx->SetTitle("");
  grPdgdEdx->SetLineWidth(2);
  grPdgdEdx->Draw("ac");
  grPdgdEdx->GetXaxis()->SetTitle("KE_{#mu} (MeV)");
  grPdgdEdx->GetYaxis()->SetTitle("MPV dE/dx (MeV/cm)");
  grVavdEdx->SetLineColor(2);
  grVavdEdx->SetLineWidth(2);
  grVavdEdx->Draw("c");
  gPad->SetLogx();
  TLegend *leg = new TLegend(0.5,0.7,0.85,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(grVavdEdx, "TMath::Vavilov","l");
  leg->AddEntry(grPdgdEdx, "PDG","l");
  leg->Draw();

  TGraph *grKappaKE = new TGraph(n2, KE, vkappa);

  TCanvas *c2 = new TCanvas("c2","c2");
  grKappaKE->SetTitle(";KE_{#mu} (MeV);#kappa");
  grKappaKE->SetLineWidth(2);
  grKappaKE->Draw("ac");
  gPad->SetLogx();
  gPad->SetLogy();

  c1->Print("MPVdEdx.pdf");
  c1->Print("MPVdEdx.png");

  c2->Print("kappaKE.pdf");
  c2->Print("kappaKE.png");
  
}
