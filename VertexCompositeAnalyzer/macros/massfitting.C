//  This is a macro for fitting mass of d0 produced at cms.
//  Also, a function drawing the invariant mass is defined, named after drawMassFitting
//
//  Here we have a short instruction for the mass fitting and drawing function
//
//  1. There are several parameters for mass fitting:
//    The 1st parameter is the histogram represneting the data               
//    The 2nd parameter is the histogram represneting the signal of MC samples 
//    The 3rd parameter is the histogram represneting 
//        the signal and swap of MC samples  
//    The 4th parameter is the name of the fitting function            
//    The 5th parameter returns the fitting result
//
//  2. Drawing function has a list of paramters:
//    1) the histogram of data
//    2) the fiting function
//    3) the name passed to TCanvas.Print(name);
//    4) whether the histogram is normalized
//    5) the bin label indicating the pT and y range
//    6) the (mva) cut applied
//    one more thing to notice is that one have to specifiy the y axis title 
//    before calling the drawing function

#ifndef __MASS_FITTING__
#define __MASS_FITTING__

#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

TF1 massfitting
(TH1* h_data, TH1* h_mc_match_signal, TH1* h_mc_match_all, const char* fname, TFitResultPtr& fitResultPtr)
{
//   The full fitting function is constructed as follow
//   [0] is signal + swap yield;
//   [1] is common mean of double gaussian;
//   [2] is signal gaussian 1 sigma;
//   [3] is signal gaussian 2 sigma;
//   [4] is fractional signal gaussian 1 yield; 
//   1-[4] is fractional signal gaussian 2 yield;
//   [5] is fractional double gaussian signal yield, 
//   1-[5] is fractional swap yield;
//   [6] is a factor to let width of the gaussians to vary in data;
//   [7] is swap gaussian sigma;
//   [8] is swap gaussian mean;
//   [9-12] is 3rd order poly parameters
    const std::string double_gaussian = 
        "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))"
        "/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
        "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))"
        "/(sqrt(2*3.14159)*[3]*(1.0 +[6]))";
    const std::string bkg3rdpoly = "[9] + [10]*x + [11]*x*x + [12]*x*x*x";
    const std::string massfunc = "[0]*(" + double_gaussian +")"
        " + " + bkg3rdpoly;
    double fit_range_low = 1.045;
    double fit_range_high = 2.0;
    double D0_mass = 1.8648;
    TF1 f(fname, massfunc.c_str(), fit_range_low, fit_range_high);

    //first fit MC signal, swap and poly bkg set to 0
        
    f.SetParameter(0,100.);
    f.SetParameter(1,D0_mass);
    f.SetParameter(2,0.03);
    f.SetParameter(3,0.005);
    f.SetParameter(4,0.1);
    
    f.FixParameter(5,1);
    f.FixParameter(6,0); //always 0 in MC
    f.FixParameter(7,0.1); //does not really mater here as yield is fix to 0
    f.FixParameter(8,D0_mass); //does not really mater here 
                                // as yield is fix to 0
    f.FixParameter(9,0);
    f.FixParameter(10,0);
    f.FixParameter(11,0);
    f.FixParameter(12,0);
    
    f.SetParLimits(2,0.01,0.1);
    f.SetParLimits(3,0.001,0.05);
    f.SetParLimits(4,0,1);
    f.SetParLimits(5,0,1);
    
    f.FixParameter(1,1.8648); //for first few attempt fix mean of gaussian 
                               //to get reasonable estimation of other pars;
                               // later open it up
    h_mc_match_signal->Fit(&f,"q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"q","",fit_range_low,fit_range_high);
    f.ReleaseParameter(1); //now let gaussian mean float
    h_mc_match_signal->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_signal->Fit(&f,"L m q","",fit_range_low,fit_range_high);

    //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
    f.FixParameter(1,f.GetParameter(1));
    f.FixParameter(2,f.GetParameter(2));
    f.FixParameter(3,f.GetParameter(3));
    f.FixParameter(4,f.GetParameter(4));
    
    //now release swap bkg parameters to fit signal+swap MC
    f.ReleaseParameter(5);
    f.ReleaseParameter(7);
    f.ReleaseParameter(8);
    
    f.SetParameter(7,0.1);
    f.SetParameter(8,D0_mass);
    
    //fit signal+swap MC
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_mc_match_all->Fit(&f,"L m q","",fit_range_low,fit_range_high);


    //now fix swap bkg parameters to fit data
    f.FixParameter(5,f.GetParameter(5));
    f.FixParameter(7,f.GetParameter(7));
    f.FixParameter(8,f.GetParameter(8));

    //now release poly bkg pars
    f.ReleaseParameter(9);
    f.ReleaseParameter(10);
    f.ReleaseParameter(11);
    f.ReleaseParameter(12);
    
    //now fit data
//    f.SetParLimits(0, 0, 1e8);
    f.SetParLimits(1, 1.86, 1.87);

    h_data->Fit(&f,"q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"q","",fit_range_low,fit_range_high);

    f.ReleaseParameter(0);
    //f.ReleaseParameter(1);
    h_data->Fit(&f,"q","",fit_range_low,fit_range_high);

    f.ReleaseParameter(1); //allow data to have different 
                            //  mass peak mean than MC
    f.ReleaseParameter(6); //allow data to have different peak width than MC
    f.SetParameter(6,0);
    //f.SetParLimits(6,-1,1);
    f.SetParLimits(6,-0.5,0.5);
    //f.FixParameter(5,1);
    //
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    h_data->Fit(&f,"L q","",fit_range_low,fit_range_high);
    fitResultPtr = h_data->Fit(&f,"L m q S","",fit_range_low,fit_range_high);

    return f;
}

void drawMassFitting(TH1* hMassData, const TF1& f, const char* picName, const bool& isNormalized,
      const char* binLabel = "", const char* cut = "")
{
   TCanvas cMass("cMass", "", 550, 450);
   cMass.SetLeftMargin(0.16);
   cMass.SetBottomMargin(0.16);

   TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
         "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
   signal.SetLineColor(kOrange-3);
   signal.SetLineWidth(1);
   signal.SetLineStyle(2);
   signal.SetFillColorAlpha(kOrange-3,0.3);
   signal.SetFillStyle(1001);
   for(int ipar=0; ipar<6+1; ipar++){
      signal.FixParameter(ipar, f.GetParameter(ipar));
   }

   TF1 swap("swap", "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))"
         "+0 *[1]*[2]*[3]*[4]", 1.7, 2.0);
   swap.SetLineColor(kGreen+4);
   swap.SetLineWidth(1);
   swap.SetLineStyle(1);
   swap.SetFillColorAlpha(kGreen+4,0.3);
   swap.SetFillStyle(1001);
   for(int ipar=0; ipar<8+1; ipar++){
      swap.FixParameter(ipar, f.GetParameter(ipar));
   }

   TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
         "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
   bkg.SetLineColor(4);
   bkg.SetLineWidth(1);
   bkg.SetLineStyle(2);
   for(int ipar=0; ipar<12+1; ipar++){
      bkg.FixParameter(ipar, f.GetParameter(ipar));
   }

   hMassData->SetTitle("");
   hMassData->GetXaxis()->SetTitle("Mass (GeV)");
   hMassData->SetTitleOffset(1.2, "y");
   hMassData->GetYaxis()->SetRangeUser(0, 1.3*hMassData->GetMaximum());
   auto gYaxis = (TGaxis*) hMassData->GetYaxis();
   gYaxis->SetMaxDigits(3);
   hMassData->Draw("E P");
   signal.Draw("same FC");
   swap.Draw("same FC");
   bkg.Draw("same L");
   hMassData->SetMarkerStyle(20);
   hMassData->Draw("SAME E P");

   TLegend lgdMass(0.6, 0.7, 0.9, 0.9);
   lgdMass.AddEntry(hMassData, "Data", "p");
   lgdMass.AddEntry(&f, "Fitting", "l");
   lgdMass.AddEntry(&signal, "Signal", "f");
   lgdMass.AddEntry(&swap, "Swap", "f");
   lgdMass.AddEntry(&bkg, "Background", "l");
   lgdMass.Draw();

   TLatex ltx;
   ltx.SetTextFont(42);
   ltx.SetTextSize(0.035);
   ltx.DrawLatexNDC(0.25, 0.8, cut);
   ltx.DrawLatexNDC(0.2, 0.65, binLabel);
   if(isNormalized) {
      ltx.DrawLatexNDC(0.65, 0.65, Form("s+swap: \n%.1f +/- %.1f", f.GetParameter(0), f.GetParError(0)));
   }else{
      ltx.DrawLatexNDC(0.65, 0.65, Form("s+swap = \n%.1f +/- %.1f", 
               f.GetParameter(0)/hMassData->GetBinWidth(1), f.GetParError(0)/hMassData->GetBinWidth(1)));
   }

   cMass.Print(picName);

}
#endif
