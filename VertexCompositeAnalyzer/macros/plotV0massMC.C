//#include "massfitting.C"

#define ks_mass 0.497611
#define lam_mass 1.115683
#define xi_mass 1.32222
#define om_mass 1.673

void plotV0massMC()
{
  TFile* f = new TFile("jetv0MC_output_08122021.root");

  TH2D* hLamRecoAll_trkdca1 = (TH2D*)f->Get("hLamRecoAll_trkdca1");
  TH2D* hLamRecoAll_trkdca2 = (TH2D*)f->Get("hLamRecoAll_trkdca2");
  TH2D* hLamRecoAll_trkdcaxy1 = (TH2D*)f->Get("hLamRecoAll_trkdcaxy1");
  TH2D* hLamRecoAll_trkdcaxy2 = (TH2D*)f->Get("hLamRecoAll_trkdcaxy2");
  TH2D* hLamRecoAll_trkdcaz1 = (TH2D*)f->Get("hLamRecoAll_trkdcaz1");
  TH2D* hLamRecoAll_trkdcaz2 = (TH2D*)f->Get("hLamRecoAll_trkdcaz2");
  TH2D* hLamRecoAll_angle3D = (TH2D*)f->Get("hLamRecoAll_angle3D");
  TH2D* hLamRecoAll_dl3D = (TH2D*)f->Get("hLamRecoAll_dl3D");
  TH2D* hLamRecoMatch_trkdca1 = (TH2D*)f->Get("hLamRecoMatch_trkdca1");
  TH2D* hLamRecoMatch_trkdca2 = (TH2D*)f->Get("hLamRecoMatch_trkdca2");
  TH2D* hLamRecoMatch_trkdcaxy1 = (TH2D*)f->Get("hLamRecoMatch_trkdcaxy1");
  TH2D* hLamRecoMatch_trkdcaxy2 = (TH2D*)f->Get("hLamRecoMatch_trkdcaxy2");
  TH2D* hLamRecoMatch_trkdcaz1 = (TH2D*)f->Get("hLamRecoMatch_trkdcaz1");
  TH2D* hLamRecoMatch_trkdcaz2 = (TH2D*)f->Get("hLamRecoMatch_trkdcaz2");
  TH2D* hLamRecoMatch_angle3D = (TH2D*)f->Get("hLamRecoMatch_angle3D");
  TH2D* hLamRecoMatch_dl3D = (TH2D*)f->Get("hLamRecoMatch_dl3D");

  TH2D* hLamRecoUnMatch_trkdca1 = (TH2D*)hLamRecoAll_trkdca1->Clone("hLamRecoUnMatch_trkdca1");
  TH2D* hLamRecoUnMatch_trkdca2 = (TH2D*)hLamRecoAll_trkdca2->Clone("hLamRecoUnMatch_trkdca2");
  TH2D* hLamRecoUnMatch_trkdcaxy1 = (TH2D*)hLamRecoAll_trkdcaxy1->Clone("hLamRecoUnMatch_trkdcaxy1");
  TH2D* hLamRecoUnMatch_trkdcaxy2 = (TH2D*)hLamRecoAll_trkdcaxy2->Clone("hLamRecoUnMatch_trkdcaxy2");
  TH2D* hLamRecoUnMatch_trkdcaz1 = (TH2D*)hLamRecoAll_trkdcaz1->Clone("hLamRecoUnMatch_trkdcaz1");
  TH2D* hLamRecoUnMatch_trkdcaz2 = (TH2D*)hLamRecoAll_trkdcaz2->Clone("hLamRecoUnMatch_trkdcaz2");
  TH2D* hLamRecoUnMatch_angle3D = (TH2D*)hLamRecoAll_angle3D->Clone("hLamRecoUnMatch_angle3D");
  TH2D* hLamRecoUnMatch_dl3D = (TH2D*)hLamRecoAll_dl3D->Clone("hLamRecoUnMatch_dl3D");

  hLamRecoUnMatch_trkdca1->Add(hLamRecoMatch_trkdca1,-1);
  hLamRecoUnMatch_trkdca2->Add(hLamRecoMatch_trkdca2,-1);
  hLamRecoUnMatch_trkdcaxy1->Add(hLamRecoMatch_trkdcaxy1,-1);
  hLamRecoUnMatch_trkdcaxy2->Add(hLamRecoMatch_trkdcaxy2,-1);
  hLamRecoUnMatch_trkdcaz1->Add(hLamRecoMatch_trkdcaz1,-1);
  hLamRecoUnMatch_trkdcaz2->Add(hLamRecoMatch_trkdcaz2,-1);
  hLamRecoUnMatch_angle3D->Add(hLamRecoMatch_angle3D,-1);
  hLamRecoUnMatch_dl3D->Add(hLamRecoMatch_dl3D,-1);

  TH2D* hXiRecoAll_trkdca1 = (TH2D*)f->Get("hXiRecoAll_trkdca1");
  TH2D* hXiRecoAll_gtrkdca1 = (TH2D*)f->Get("hXiRecoAll_gtrkdca1");
  TH2D* hXiRecoAll_gtrkdca2 = (TH2D*)f->Get("hXiRecoAll_gtrkdca2");
  TH2D* hXiRecoAll_lamangle3D = (TH2D*)f->Get("hXiRecoAll_lamangle3D");
  TH2D* hXiRecoAll_lamdl3D = (TH2D*)f->Get("hXiRecoAll_lamdl3D");
  TH2D* hXiRecoAll_angle3D = (TH2D*)f->Get("hXiRecoAll_angle3D");
  TH2D* hXiRecoAll_dl3D = (TH2D*)f->Get("hXiRecoAll_dl3D");
  TH2D* hXiRecoMatch_trkdca1 = (TH2D*)f->Get("hXiRecoMatch_trkdca1");
  TH2D* hXiRecoMatch_gtrkdca1 = (TH2D*)f->Get("hXiRecoMatch_gtrkdca1");
  TH2D* hXiRecoMatch_gtrkdca2 = (TH2D*)f->Get("hXiRecoMatch_gtrkdca2");
  TH2D* hXiRecoMatch_lamangle3D = (TH2D*)f->Get("hXiRecoMatch_lamangle3D");
  TH2D* hXiRecoMatch_lamdl3D = (TH2D*)f->Get("hXiRecoMatch_lamdl3D");
  TH2D* hXiRecoMatch_angle3D = (TH2D*)f->Get("hXiRecoMatch_angle3D");
  TH2D* hXiRecoMatch_dl3D = (TH2D*)f->Get("hXiRecoMatch_dl3D");

  TH2D* hXiRecoUnMatch_trkdca1 = (TH2D*)hXiRecoAll_trkdca1->Clone("hXiRecoUnMatch_trkdca1");
  TH2D* hXiRecoUnMatch_gtrkdca1 = (TH2D*)hXiRecoAll_gtrkdca1->Clone("hXiRecoUnMatch_gtrkdca1");
  TH2D* hXiRecoUnMatch_gtrkdca2 = (TH2D*)hXiRecoAll_gtrkdca2->Clone("hXiRecoUnMatch_gtrkdca2");
  TH2D* hXiRecoUnMatch_lamangle3D = (TH2D*)hXiRecoAll_lamangle3D->Clone("hXiRecoUnMatch_lamangle3D");
  TH2D* hXiRecoUnMatch_lamdl3D = (TH2D*)hXiRecoAll_lamdl3D->Clone("hXiRecoUnMatch_lamdl3D");
  TH2D* hXiRecoUnMatch_angle3D = (TH2D*)hXiRecoAll_angle3D->Clone("hXiRecoUnMatch_angle3D");
  TH2D* hXiRecoUnMatch_dl3D = (TH2D*)hXiRecoAll_dl3D->Clone("hXiRecoUnMatch_dl3D");

  hXiRecoUnMatch_trkdca1->Add(hXiRecoMatch_trkdca1,-1);
  hXiRecoUnMatch_gtrkdca1->Add(hXiRecoMatch_gtrkdca1,-1);
  hXiRecoUnMatch_gtrkdca2->Add(hXiRecoMatch_gtrkdca2,-1);
  hXiRecoUnMatch_lamangle3D->Add(hXiRecoMatch_lamangle3D,-1);
  hXiRecoUnMatch_lamdl3D->Add(hXiRecoMatch_lamdl3D,-1);
  hXiRecoUnMatch_angle3D->Add(hXiRecoMatch_angle3D,-1);
  hXiRecoUnMatch_dl3D->Add(hXiRecoMatch_dl3D,-1);

  int istart=2;
  int iend=10;
  TH1D* hLamRecoUnMatch1D_trkdca1 = (TH1D*)hLamRecoUnMatch_trkdca1->ProjectionX("hLamRecoUnMatch1D_trkdca1",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_trkdca2 = (TH1D*)hLamRecoUnMatch_trkdca2->ProjectionX("hLamRecoUnMatch1D_trkdca2",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_trkdcaxy1 = (TH1D*)hLamRecoUnMatch_trkdcaxy1->ProjectionX("hLamRecoUnMatch1D_trkdcaxy1",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_trkdcaxy2 = (TH1D*)hLamRecoUnMatch_trkdcaxy2->ProjectionX("hLamRecoUnMatch1D_trkdcaxy2",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_trkdcaz1 = (TH1D*)hLamRecoUnMatch_trkdcaz1->ProjectionX("hLamRecoUnMatch1D_trkdcaz1",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_trkdcaz2 = (TH1D*)hLamRecoUnMatch_trkdcaz2->ProjectionX("hLamRecoUnMatch1D_trkdcaz2",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_angle3D = (TH1D*)hLamRecoUnMatch_angle3D->ProjectionX("hLamRecoUnMatch1D_angle3D",istart,iend,"e");
  TH1D* hLamRecoUnMatch1D_dl3D = (TH1D*)hLamRecoUnMatch_dl3D->ProjectionX("hLamRecoUnMatch1D_dl3D",istart,iend,"e");

  TH1D* hLamRecoMatch1D_trkdca1 = (TH1D*)hLamRecoMatch_trkdca1->ProjectionX("hLamRecoMatch1D_trkdca1",istart,iend,"e");
  TH1D* hLamRecoMatch1D_trkdca2 = (TH1D*)hLamRecoMatch_trkdca2->ProjectionX("hLamRecoMatch1D_trkdca2",istart,iend,"e");
  TH1D* hLamRecoMatch1D_trkdcaxy1 = (TH1D*)hLamRecoMatch_trkdcaxy1->ProjectionX("hLamRecoMatch1D_trkdcaxy1",istart,iend,"e");
  TH1D* hLamRecoMatch1D_trkdcaxy2 = (TH1D*)hLamRecoMatch_trkdcaxy2->ProjectionX("hLamRecoMatch1D_trkdcaxy2",istart,iend,"e");
  TH1D* hLamRecoMatch1D_trkdcaz1 = (TH1D*)hLamRecoMatch_trkdcaz1->ProjectionX("hLamRecoMatch1D_trkdcaz1",istart,iend,"e");
  TH1D* hLamRecoMatch1D_trkdcaz2 = (TH1D*)hLamRecoMatch_trkdcaz2->ProjectionX("hLamRecoMatch1D_trkdcaz2",istart,iend,"e");
  TH1D* hLamRecoMatch1D_angle3D = (TH1D*)hLamRecoMatch_angle3D->ProjectionX("hLamRecoMatch1D_angle3D",istart,iend,"e");
  TH1D* hLamRecoMatch1D_dl3D = (TH1D*)hLamRecoMatch_dl3D->ProjectionX("hLamRecoMatch1D_dl3D",istart,iend,"e");

  hLamRecoUnMatch1D_trkdca1->SetMarkerStyle(24);
  hLamRecoUnMatch1D_trkdca2->SetMarkerStyle(24);
  hLamRecoUnMatch1D_trkdcaxy1->SetMarkerStyle(24);
  hLamRecoUnMatch1D_trkdcaxy2->SetMarkerStyle(24);
  hLamRecoUnMatch1D_trkdcaz1->SetMarkerStyle(24);
  hLamRecoUnMatch1D_trkdcaz2->SetMarkerStyle(24);
  hLamRecoUnMatch1D_angle3D->SetMarkerStyle(24);
  hLamRecoUnMatch1D_dl3D->SetMarkerStyle(24);

  hLamRecoMatch1D_trkdca1->SetMarkerStyle(20);
  hLamRecoMatch1D_trkdca2->SetMarkerStyle(20);
  hLamRecoMatch1D_trkdcaxy1->SetMarkerStyle(20);
  hLamRecoMatch1D_trkdcaxy2->SetMarkerStyle(20);
  hLamRecoMatch1D_trkdcaz1->SetMarkerStyle(20);
  hLamRecoMatch1D_trkdcaz2->SetMarkerStyle(20);
  hLamRecoMatch1D_angle3D->SetMarkerStyle(20);
  hLamRecoMatch1D_dl3D->SetMarkerStyle(20);

  hLamRecoUnMatch1D_trkdca1->Scale(1.0/hLamRecoUnMatch1D_trkdca1->Integral("width"));
  hLamRecoUnMatch1D_trkdca2->Scale(1.0/hLamRecoUnMatch1D_trkdca2->Integral("width"));
  hLamRecoUnMatch1D_trkdcaxy1->Scale(1.0/hLamRecoUnMatch1D_trkdcaxy1->Integral("width"));
  hLamRecoUnMatch1D_trkdcaxy2->Scale(1.0/hLamRecoUnMatch1D_trkdcaxy2->Integral("width"));
  hLamRecoUnMatch1D_trkdcaz1->Scale(1.0/hLamRecoUnMatch1D_trkdcaz1->Integral("width"));
  hLamRecoUnMatch1D_trkdcaz2->Scale(1.0/hLamRecoUnMatch1D_trkdcaz2->Integral("width"));
  hLamRecoUnMatch1D_angle3D->Scale(1.0/hLamRecoUnMatch1D_angle3D->Integral("width"));
  hLamRecoUnMatch1D_dl3D->Scale(1.0/hLamRecoUnMatch1D_dl3D->Integral("width"));

  hLamRecoMatch1D_trkdca1->Rebin(5);
  hLamRecoMatch1D_trkdca2->Rebin(5);
  hLamRecoMatch1D_trkdcaxy1->Rebin(5);
  hLamRecoMatch1D_trkdcaxy2->Rebin(5);
  hLamRecoMatch1D_trkdcaz1->Rebin(5);
  hLamRecoMatch1D_trkdcaz2->Rebin(5);
  hLamRecoMatch1D_dl3D->Rebin(5);

  hLamRecoMatch1D_trkdca1->Scale(1.0/hLamRecoMatch1D_trkdca1->Integral("width"));
  hLamRecoMatch1D_trkdca2->Scale(1.0/hLamRecoMatch1D_trkdca2->Integral("width"));
  hLamRecoMatch1D_trkdcaxy1->Scale(1.0/hLamRecoMatch1D_trkdcaxy1->Integral("width"));
  hLamRecoMatch1D_trkdcaxy2->Scale(1.0/hLamRecoMatch1D_trkdcaxy2->Integral("width"));
  hLamRecoMatch1D_trkdcaz1->Scale(1.0/hLamRecoMatch1D_trkdcaz1->Integral("width"));
  hLamRecoMatch1D_trkdcaz2->Scale(1.0/hLamRecoMatch1D_trkdcaz2->Integral("width"));
  hLamRecoMatch1D_angle3D->Scale(1.0/hLamRecoMatch1D_angle3D->Integral("width"));
  hLamRecoMatch1D_dl3D->Scale(1.0/hLamRecoMatch1D_dl3D->Integral("width"));

  TH1D* hXiRecoUnMatch1D_trkdca1 = (TH1D*)hXiRecoUnMatch_trkdca1->ProjectionX("hXiRecoUnMatch1D_trkdca1",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_gtrkdca1 = (TH1D*)hXiRecoUnMatch_gtrkdca1->ProjectionX("hXiRecoUnMatch1D_gtrkdca1",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_gtrkdca2 = (TH1D*)hXiRecoUnMatch_gtrkdca2->ProjectionX("hXiRecoUnMatch1D_gtrkdca2",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_lamangle3D = (TH1D*)hXiRecoUnMatch_lamangle3D->ProjectionX("hXiRecoUnMatch1D_lamangle3D",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_lamdl3D = (TH1D*)hXiRecoUnMatch_lamdl3D->ProjectionX("hXiRecoUnMatch1D_lamdl3D",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_angle3D = (TH1D*)hXiRecoUnMatch_angle3D->ProjectionX("hXiRecoUnMatch1D_angle3D",istart,iend,"e");
  TH1D* hXiRecoUnMatch1D_dl3D = (TH1D*)hXiRecoUnMatch_dl3D->ProjectionX("hXiRecoUnMatch1D_dl3D",istart,iend,"e");

  TH1D* hXiRecoMatch1D_trkdca1 = (TH1D*)hXiRecoMatch_trkdca1->ProjectionX("hXiRecoMatch1D_trkdca1",istart,iend,"e");
  TH1D* hXiRecoMatch1D_gtrkdca1 = (TH1D*)hXiRecoMatch_gtrkdca1->ProjectionX("hXiRecoMatch1D_gtrkdca1",istart,iend,"e");
  TH1D* hXiRecoMatch1D_gtrkdca2 = (TH1D*)hXiRecoMatch_gtrkdca2->ProjectionX("hXiRecoMatch1D_gtrkdca2",istart,iend,"e");
  TH1D* hXiRecoMatch1D_lamangle3D = (TH1D*)hXiRecoMatch_lamangle3D->ProjectionX("hXiRecoMatch1D_lamangle3D",istart,iend,"e");
  TH1D* hXiRecoMatch1D_lamdl3D = (TH1D*)hXiRecoMatch_lamdl3D->ProjectionX("hXiRecoMatch1D_lamdl3D",istart,iend,"e");
  TH1D* hXiRecoMatch1D_angle3D = (TH1D*)hXiRecoMatch_angle3D->ProjectionX("hXiRecoMatch1D_angle3D",istart,iend,"e");
  TH1D* hXiRecoMatch1D_dl3D = (TH1D*)hXiRecoMatch_dl3D->ProjectionX("hXiRecoMatch1D_dl3D",istart,iend,"e");

  hXiRecoUnMatch1D_trkdca1->SetMarkerStyle(24);
  hXiRecoUnMatch1D_gtrkdca1->SetMarkerStyle(24);
  hXiRecoUnMatch1D_gtrkdca2->SetMarkerStyle(24);
  hXiRecoUnMatch1D_lamangle3D->SetMarkerStyle(24);
  hXiRecoUnMatch1D_lamdl3D->SetMarkerStyle(24);
  hXiRecoUnMatch1D_angle3D->SetMarkerStyle(24);
  hXiRecoUnMatch1D_dl3D->SetMarkerStyle(24);

  hXiRecoMatch1D_trkdca1->SetMarkerStyle(20);
  hXiRecoMatch1D_gtrkdca1->SetMarkerStyle(20);
  hXiRecoMatch1D_gtrkdca2->SetMarkerStyle(20);
  hXiRecoMatch1D_lamangle3D->SetMarkerStyle(20);
  hXiRecoMatch1D_lamdl3D->SetMarkerStyle(20);
  hXiRecoMatch1D_angle3D->SetMarkerStyle(20);
  hXiRecoMatch1D_dl3D->SetMarkerStyle(20);

  hXiRecoUnMatch1D_trkdca1->Rebin(10);
  hXiRecoUnMatch1D_gtrkdca1->Rebin(10);
  hXiRecoUnMatch1D_gtrkdca2->Rebin(10);
  hXiRecoUnMatch1D_lamdl3D->Rebin(10);
  hXiRecoUnMatch1D_dl3D->Rebin(10);

  hXiRecoUnMatch1D_trkdca1->Scale(1.0/hXiRecoUnMatch1D_trkdca1->Integral("width"));
  hXiRecoUnMatch1D_gtrkdca1->Scale(1.0/hXiRecoUnMatch1D_gtrkdca1->Integral("width"));
  hXiRecoUnMatch1D_gtrkdca2->Scale(1.0/hXiRecoUnMatch1D_gtrkdca2->Integral("width"));
  hXiRecoUnMatch1D_lamangle3D->Scale(1.0/hXiRecoUnMatch1D_lamangle3D->Integral("width"));
  hXiRecoUnMatch1D_lamdl3D->Scale(1.0/hXiRecoUnMatch1D_lamdl3D->Integral("width"));
  hXiRecoUnMatch1D_angle3D->Scale(1.0/hXiRecoUnMatch1D_angle3D->Integral("width"));
  hXiRecoUnMatch1D_dl3D->Scale(1.0/hXiRecoUnMatch1D_dl3D->Integral("width"));

  hXiRecoMatch1D_trkdca1->Rebin(20);
  hXiRecoMatch1D_gtrkdca1->Rebin(20);
  hXiRecoMatch1D_gtrkdca2->Rebin(20);
  hXiRecoMatch1D_lamdl3D->Rebin(20);
  hXiRecoMatch1D_dl3D->Rebin(20);

  hXiRecoMatch1D_trkdca1->Scale(1.0/hXiRecoMatch1D_trkdca1->Integral("width"));
  hXiRecoMatch1D_gtrkdca1->Scale(1.0/hXiRecoMatch1D_gtrkdca1->Integral("width"));
  hXiRecoMatch1D_gtrkdca2->Scale(1.0/hXiRecoMatch1D_gtrkdca2->Integral("width"));
  hXiRecoMatch1D_lamangle3D->Scale(1.0/hXiRecoMatch1D_lamangle3D->Integral("width"));
  hXiRecoMatch1D_lamdl3D->Scale(1.0/hXiRecoMatch1D_lamdl3D->Integral("width"));
  hXiRecoMatch1D_angle3D->Scale(1.0/hXiRecoMatch1D_angle3D->Integral("width"));
  hXiRecoMatch1D_dl3D->Scale(1.0/hXiRecoMatch1D_dl3D->Integral("width"));


  TCanvas* ccc = new TCanvas("ccc","ccc",1600,900);
  ccc->Divide(4,2);
  ccc->cd(1);
  hLamRecoMatch1D_trkdca1->Draw("PE");
  hLamRecoUnMatch1D_trkdca1->Draw("PESAME");
  ccc->cd(2);
  hLamRecoMatch1D_trkdca2->Draw("PE");
  hLamRecoUnMatch1D_trkdca2->Draw("PESAME");
  ccc->cd(3);
  hLamRecoMatch1D_trkdcaxy1->Draw("PE");
  hLamRecoUnMatch1D_trkdcaxy1->Draw("PESAME");
  ccc->cd(4);
  hLamRecoMatch1D_trkdcaxy2->Draw("PE");
  hLamRecoUnMatch1D_trkdcaxy2->Draw("PESAME");
  ccc->cd(5);
  hLamRecoMatch1D_trkdcaz1->Draw("PE");
  hLamRecoUnMatch1D_trkdcaz1->Draw("PESAME");
  ccc->cd(6);
  hLamRecoMatch1D_trkdcaz2->Draw("PE");
  hLamRecoUnMatch1D_trkdcaz2->Draw("PESAME");
  ccc->cd(7);
  hLamRecoMatch1D_angle3D->SetAxisRange(0.999,1,"X");
  hLamRecoMatch1D_angle3D->Draw("PE");
  hLamRecoUnMatch1D_angle3D->Draw("PESAME");
  ccc->cd(8);
  hLamRecoMatch1D_dl3D->Draw("PE");
  hLamRecoUnMatch1D_dl3D->Draw("PESAME");


  TCanvas* ccc1 = new TCanvas("ccc1","ccc1",1600,900);
  ccc1->Divide(4,2);
  ccc1->cd(1);
  hXiRecoMatch1D_trkdca1->Draw("PE");
  hXiRecoUnMatch1D_trkdca1->Draw("PESAME");
  ccc1->cd(2);
  hXiRecoMatch1D_gtrkdca1->Draw("PE");
  hXiRecoUnMatch1D_gtrkdca1->Draw("PESAME");
  ccc1->cd(3);
  hXiRecoMatch1D_gtrkdca2->Draw("PE");
  hXiRecoUnMatch1D_gtrkdca2->Draw("PESAME");
  ccc1->cd(4);
  hXiRecoMatch1D_lamangle3D->SetAxisRange(0.999,1,"X");
  hXiRecoMatch1D_lamangle3D->Draw("PE");
  hXiRecoUnMatch1D_lamangle3D->Draw("PESAME");
  ccc1->cd(5);
  hXiRecoMatch1D_lamdl3D->Draw("PE");
  hXiRecoUnMatch1D_lamdl3D->Draw("PESAME");
  ccc1->cd(6);
  hXiRecoMatch1D_angle3D->SetAxisRange(0.999,1,"X");
  hXiRecoMatch1D_angle3D->Draw("PE");
  hXiRecoUnMatch1D_angle3D->Draw("PESAME");
  ccc1->cd(7);
  hXiRecoMatch1D_dl3D->Draw("PE");
  hXiRecoUnMatch1D_dl3D->Draw("PESAME");

return;
  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8OmRecoAll_mass");
//  TH2D* hAntiLamMassVsPt = (TH2D*)f->Get("hJet8AntiXiRecoAll_mass");
//  hLamMassVsPt->Add(hAntiLamMassVsPt);
//  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8XiRecoAll_mass");
//  TH2D* hLamMassVsPt = (TH2D*)f->Get("hJet8KsRecoAll_mass");
  TH1D* hLamMass[20];

//  for(int i=0;i<hLamMassVsPt->GetNbinsX();i++)
  for(int i=0;i<20;i++)
//    hLamMass[i] = (TH1D*)hLamMassVsPt->ProjectionY(Form("hJetLamMass_%d",i),i+1,i+1,"e");
    hLamMass[i] = (TH1D*)hLamMassVsPt->ProjectionY(Form("hJetLamMass_%d",i),8*i+1,8*i+8,"e");

  const std::string double_gaussian =
      "[4]*TMath::Gaus(x,[1],[2])"
      "/(sqrt(2*3.14159)*[2])"
      "+ (1-[4])*TMath::Gaus(x,[1],[3])"
      "/(sqrt(2*3.14159)*[3])";
  const std::string bkg3rdpoly = "[5] + [6]*x + [7]*x*x + [8]*x*x*x";
//  const std::string bkg3rdpoly = "[5]*pow(x-0.13957-0.93827,0.5) + [6]*pow(x-0.13957-0.93827,1.5) + [7]*pow(x-0.13957-0.93827,2.5)";
//  const std::string bkg3rdpoly = "[5]*pow(x-0.13957-1.115683,0.5) + [6]*pow(x-0.13957-1.115683,1.5) + [7]*pow(x-0.13957-1.115683,2.5)";
  const std::string massfunc = "[0]*(" + double_gaussian +")"
      " + " + bkg3rdpoly;

  TF1* fitfunc[20];
  for(int i=0;i<20;i++)
  {
//   fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),0.44,0.54);
//   fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.09,1.144);
//    fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.26,1.38);
    fitfunc[i] = new TF1(Form("fitfunc_%d",i),massfunc.c_str(),1.6,1.74);
    fitfunc[i]->SetParameter(0,500);
//    fitfunc[i]->SetParameter(1,ks_mass);
//    fitfunc[i]->SetParameter(1,lam_mass); 
//    fitfunc[i]->SetParameter(1,xi_mass);
    fitfunc[i]->SetParameter(1,om_mass);
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
//    hLamMass[i]->SetAxisRange(1.27,1.37,"X");
    hLamMass[i]->Draw("PE");
    fitfunc[i]->Draw("Lsame");
  }

}
