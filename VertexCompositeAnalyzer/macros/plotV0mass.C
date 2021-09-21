//#include "massfitting.C"

#define ks_mass 0.497611
#define lam_mass 1.115683
#define xi_mass 1.32222
#define om_mass 1.673

void plotV0mass()
{
  TFile* f = new TFile("jetv0_output_set19_08072021.root");
//  TFile* f = new TFile("jetv0_output.root");
//  TFile* f = new TFile("jetv0_output_filelist_v0.txt.root");
//  TFile* f = new TFile("jetv0_output_list_v0_090.root");

  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8XiRecoAll_mass");
//  TH2D* hAntiLamMassVsPt = (TH2D*)f->Get("hJet8AntiXiRecoAll_mass");
//  hLamMassVsPt->Add(hAntiLamMassVsPt);
//  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8XiRecoAll_mass");
//  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8KsRecoAll_mass");
  TH1D* hLamMass[20];

//  for(int i=0;i<hLamMassVsPt->GetNbinsX();i++)
  for(int i=0;i<20;i++)
//    hLamMass[i] = (TH1D*)hLamMassVsPt->ProjectionY(Form("hJetLamMass_%d",i),i+1,i+1,"e");
    hLamMass[i] = (TH1D*)hLamMassVsPt->ProjectionY(Form("hJetLamMass_%d",i),4*i+1,4*i+4,"e");

  const std::string double_gaussian =
      "[4]*TMath::Gaus(x,[1],[2])"
      "/(sqrt(2*3.14159)*[2])"
      "+ (1-[4])*TMath::Gaus(x,[1],[3])"
      "/(sqrt(2*3.14159)*[3])";
//  const std::string bkg3rdpoly = "[5] + [6]*x + [7]*x*x + [8]*x*x*x";
  const std::string bkg3rdpoly = "[5]*pow(x-0.13957-0.93827,0.5) + [6]*pow(x-0.13957-0.93827,1.5) + [7]*pow(x-0.13957-0.93827,2.5)";
//  const std::string bkg3rdpoly = "[5]*pow(x-0.13957-1.115683,0.5) + [6]*pow(x-0.13957-1.115683,1.5) + [7]*pow(x-0.13957-1.115683,2.5)";
  const std::string massfunc = "[0]*(" + double_gaussian +")"
      " + " + bkg3rdpoly;

  TF1* fitfunc[20];
  for(int i=0;i<20;i++)
  {
//   fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),0.44,0.54);
   fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.09,1.144);
//    fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.26,1.38);
//    fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.6,1.74);
    fitfunc[i]->SetParameter(0,500);
//    fitfunc[i]->SetParameter(1,ks_mass);
    fitfunc[i]->SetParameter(1,lam_mass); 
//    fitfunc[i]->SetParameter(1,xi_mass);
//    fitfunc[i]->SetParameter(1,om_mass);
    fitfunc[i]->SetParameter(2,0.005);
    fitfunc[i]->SetParameter(3,0.015);  
    fitfunc[i]->SetParameter(4,0.6);
    fitfunc[i]->SetParameter(5,200);
    fitfunc[i]->SetParameter(6,10);
    fitfunc[i]->FixParameter(7,0);
    fitfunc[i]->SetParLimits(4,0,1);
//    fitfunc[i]->SetParameter(8,0);

    hLamMass[i]->Fit(Form("fitfunc_%d",i),"RNO");
    hLamMass[i]->Fit(Form("fitfunc_%d",i),"RNO");
    hLamMass[i]->Fit(Form("fitfunc_%d",i),"RNO");
    hLamMass[i]->Fit(Form("fitfunc_%d",i),"RNO");
    hLamMass[i]->Fit(Form("fitfunc_%d",i),"RNO");
  }

  TCanvas* c = new TCanvas("c","",1600,1000);
  c->Divide(5,4);
  for(int i=0;i<20;i++)
  {
    c->cd(i+1);
    hLamMass[i]->SetTitle(Form("%.2f#pm%.2f, Sig=%.1f",fitfunc[i]->GetParameter(0),fitfunc[i]->GetParError(0),fitfunc[i]->GetParameter(0)/fitfunc[i]->GetParError(0)));
//    hLamMass[i]->SetAxisRange(1.085,1.144,"X");
    hLamMass[i]->SetAxisRange(1.27,1.37,"X");
    hLamMass[i]->Draw("PE");
//    fitfunc[i]->Draw("Lsame");
  }

}
