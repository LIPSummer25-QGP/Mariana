#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip>
#include "ACCSEL.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>  
#include <RooProdPdf.h>
#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooCurve.h>
using namespace RooFit;


// Aux: take the y-value from the drawn curve at a given mass (for vertical line heights)
double getYatMass(RooPlot* frame, double mass) {
    if (auto* curve = dynamic_cast<RooCurve*>(frame->findObject("global"))) {
        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();
        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    // fallback: scan any curve (your original)
    for (int i = 0; i < frame->numItems(); ++i) {
        RooCurve* curve = dynamic_cast<RooCurve*>(frame->getObject(i));
        if (!curve) continue;
        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();
        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    return 0.0;
}

// B+ Particle with Error Function
void sbs_var() {
    const int nbins_plot = 100; // Number of bins for the plot

    double min_signal = 5.165287364;//using 3 sigma effective
    double max_signal = 5.391576360;//using 3 sigma effective
    
    double mc_sigma1 = 0.03702;
    double mc_sigma2 = 0.01609;
    double mc_c1 = 0.3358;
    //double mc_c2 = 5.5593; // Not used in this fit, but defined for completeness
    double chi2_fix=1.51;
    double xlow = 5.0;
    double xhigh = 5.8;
    double Nbkg_fix = 213124.8; // Fixed background yield from previous fit
    double Nsig_fix=40781.9;
    double Nerfc_fix=8870.7;
    double lambda_fix=-0.2134;
    double cs_fix = 1.04129; // Fixed scale constant for Gaussian resolution
    double csf_fix = 5.13460; // Fixed shifting constant
    double csc_fix = 0.04354; // Fixed scaling constant for ERFC background
    double fix_mean = 5.27843;
    double xsi = 5.177858768;   //using 4 sigma     //5.077288754; // Left sideband edge 8 sigma efective
    double xsf = 5.379001232;   //using 4 sigma     //5.47957464; // Right sideband edge 8 sigma efective
  

    double bin_width_plot = (xhigh - xlow) / nbins_plot;

    // Load ROOT file and TTree
    TFile *file = TFile::Open("DATA_MC_pp_Bmass_Buforfit.root");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Tree 'tree' not found in file." << std::endl;
        return;
    }

      // Open and process the MC file for signal region yield
    TFile *file_mc = TFile::Open("/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root");
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }
    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntKp", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        return;
    }

    // Define the mass variable and dataset
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);

   // Define ranges
    B_mass.setRange("signal_Region", xsi, xsf);
    B_mass.setRange("sidebandRegion", xsf , xhigh);
    B_mass.setRange("leftsideband", xlow, xsi);
 
    // Create unbinned RooDataSet from TTree
    B_mass.setRange("gaussRange", min_signal, max_signal);

    RooArgSet vars(B_mass);
    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, vars);
    std::cout << "Number of entries in dataset: " << dataset.numEntries() << std::endl;

    // Signal model: Double Gaussian
    RooRealVar mean("mean", "Mean", fix_mean);
    mean.setConstant(kTRUE); // Fix the mean to known value

    // MC-derived widths (FIXED constants — put your values here)
    RooRealVar sigma1_mc("sigma1_mc", "MC Sigma1", mc_sigma1); 
    sigma1_mc.setConstant(kTRUE);

    RooRealVar sigma2_mc("sigma2_mc", "MC Sigma2", mc_sigma2); 
    sigma2_mc.setConstant(kTRUE);

    RooRealVar c1("c1", "Fraction of Gaussian1", mc_c1);
    c1.setConstant(kTRUE);

    // Common positive scale (fit parameter)
    RooRealVar Cs("Cs", "Resolution scale", cs_fix);
    Cs.setConstant(kTRUE); // Fix the scale to known value

    // Effective widths = Cs * sigma_mc
    RooProduct sigma1_eff("sigma1_eff", "sigma1_eff", RooArgList(sigma1_mc, Cs));
    RooProduct sigma2_eff("sigma2_eff", "sigma2_eff", RooArgList(sigma2_mc, Cs));

    // Mixture fraction and Gaussians
    RooGaussian gauss1("gauss1", "Gaussian 1", B_mass, mean, sigma1_eff);
    RooGaussian gauss2("gauss2", "Gaussian 2", B_mass, mean, sigma2_eff);
    RooAddPdf signal("signal", "Double Gaussian Model", RooArgList(gauss1, gauss2), RooArgList(c1));
    
    
    RooRealVar Nsig("Nsig", "Signal Yield", Nsig_fix);
    Nsig.setConstant(kTRUE);

    // Background model
    RooRealVar lambda("lambda", "Lambda", lambda_fix);
    lambda.setConstant(kTRUE); // Fix the lambda to known value
    RooExponential expo("expo", "Exponential Background", B_mass, lambda);    
    
    RooRealVar Nbkg("Nbkg", "Exponential Background Yield", Nbkg_fix);
    Nbkg.setConstant(kTRUE);

    // ERFC background (for left sideband)
    RooRealVar csf("csf", "Shifting Constant", csf_fix);
    csf.setConstant(kTRUE); // Fix the constant to known value
    RooRealVar csc("csc", "Scaling Constant", csc_fix);
    csc.setConstant(kTRUE); // Fix the scaling constant to known value

    // integral form implemented via erf: 1 - erf(x)
   RooGenericPdf erfc_bkg("erfc_bkg", "1 - TMath::Erf((B_mass - csf)/csc)", RooArgList(B_mass, csf, csc));

    RooRealVar Nerfc("Nerfc", "ERFC Background Yield", Nerfc_fix);
    Nerfc.setConstant(kTRUE);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, expo, erfc_bkg), RooArgList(Nsig, Nbkg, Nerfc));



// Calculate background contribution
double frac_expo_signal = expo.createIntegral(B_mass, NormSet(B_mass), Range("signal_Region"))->getVal();
double frac_expo_sideband = expo.createIntegral(B_mass, NormSet(B_mass), Range("rightsideband"))->getVal();
double bkg_signal = Nbkg.getVal() * frac_expo_signal;// bkg signal region
double bkg_sideband = Nbkg.getVal() * frac_expo_sideband;//bkg sideband region

// Number of bins for all histograms
int nbins_var = 60;
int nbins_mc = 60;
int nbins = 100;//bins for zoom plot
// Variables and their ranges
//const char *variables_data[] = {"B_mass","B_alpha","B_cos_dtheta","B_chi2cl","B_trk1dR","B_trk1Pt","B_norm_svpvDistance_2D","B_norm_svpvDistance","B_norm_trk1Dxy","B_pt","B_y"/*,"nSelectedChargedTracks"*/};

//const char *variables_mc[] = {"Bmass","Balpha","Bcos_dtheta","Bchi2cl","Btrk1dR","Btrk1Pt","Bnorm_svpvDistance_2D","Bnorm_svpvDistance","Bnorm_trk1Dxy","Bpt","By"/*,"nSelectedChargedTracks"*/};

//const double ranges[][2] = {{5,6},{0,3.15},{0,1},{0,1},{0,4.5},{0.5,10},{0,85},{0,85},{-22,22},{5,50},{-2.4,2.4}/*,{0,150}*/};



const char *variables_data[] = {"B_mass","B_alpha","B_cos_dtheta","B_trk1dR","B_trk1Pt","B_norm_svpvDistance_2D","B_norm_svpvDistance"/*,"nSelectedChargedTracks"*/};
const char *variables_mc[] = {"Bmass","Balpha","Bcos_dtheta","Btrk1dR","Btrk1Pt","Bnorm_svpvDistance_2D","Bnorm_svpvDistance"/*,"nSelectedChargedTracks"*/};

const double ranges[][2] = {{5,6},{0,0.2},{0.98,1},{1.5,4.5},{0.5,2},{0,15},{0,15}/*,{0,150}*/};





    TString cut_mc = Form(" Bnorm_svpvDistance>2 && Bnorm_svpvDistance_2D>4 && Bchi2cl>0.003 && (%s) && (%s) && (%s) && (%s)",//Btrk1dR<1.66248&& Balpha<0.150255
                        isMCsignal.Data(),
                        ACCcuts_ppRef_Bu.Data(),
                        SELcuts_ppRef_Bu.Data(),
                        TRGmatching.Data());


const int nVars = sizeof(variables_data)/sizeof(variables_data[0]);
const int nVars_mc = sizeof(variables_mc)/sizeof(variables_mc[0]);
// Sideband subtraction factor
double alpha = (bkg_sideband > 0) ? (bkg_signal / bkg_sideband) : 0.0;
std::cout << "Alpha (scaling factor): " << alpha << std::endl;
std::cout << "data tree: "  << std::endl;

for (int i = 0; i < nVars; ++i) {
    std::string var = variables_data[i];
    std::string var_mc = variables_mc[i];
    double xmin = ranges[i][0];
    double xmax = ranges[i][1];

    // Skip B_mass since that's your fit variable
    if (var == "B_mass") continue;
    if (var_mc == "Bmass") continue;

    // Create histograms for signal and sideband regions
    TH1F* hist_signal = new TH1F(Form("hist_signal_%s", var.c_str()),Form("B^{+}; %s; Entries", var.c_str()),nbins_var, xmin, xmax);
    TH1F* hist_sideband = new TH1F(Form("hist_sideband_%s", var.c_str()),Form("B^{+}; %s; Entries", var.c_str()),nbins_var, xmin, xmax);

    // Fill histograms from TTree
    tree->Draw(Form("%s >> hist_signal_%s", var.c_str(), var.c_str()),Form("B_mass >= %f && B_mass <= %f", xsi, xsf), "goff");
    tree->Draw(Form("%s >> hist_sideband_%s", var.c_str(), var.c_str()),Form("B_mass < %f || B_mass > %f", xsi, xsf), "goff");

    double Nsig = hist_signal->GetEntries();
    double Nsideband = hist_sideband->GetEntries();

    TH1F* hist_signal_subtracted = (TH1F*)hist_signal->Clone(Form("hist_signal_subtracted_%s", var.c_str()));
        hist_signal_subtracted->Add(hist_sideband, -alpha);
    
    double Nsubtracted = hist_signal_subtracted->GetEntries();

    TH1F* hist_mc = new TH1F(Form("hist_mc_%s", var_mc.c_str()),Form("%s in MC; %s; Entries", var_mc.c_str(), var_mc.c_str()),nbins_mc, xmin, xmax);
    treemc->Draw(Form("%s >> hist_mc_%s", var_mc.c_str(), var_mc.c_str()), cut_mc.Data(), "goff");
    

    if(hist_signal_subtracted->Integral()>0){ hist_signal_subtracted->Scale(1.0 / hist_signal_subtracted->Integral());}
    if (hist_mc->Integral() > 0){hist_mc->Scale(1.0 / hist_mc->Integral());}

    //Plot histograms
    TCanvas* c_compare = new TCanvas(Form("c_%s", var_mc.c_str()), var_mc.c_str(), 800, 600);  // sideband subtraction plots 
    c_compare->SetGrid();    
    c_compare->SetTicks();
        hist_mc->SetStats(0);
        hist_signal_subtracted->SetStats(0);
        hist_mc->SetMarkerStyle(20);
        hist_mc->SetMarkerSize(0.9);
        hist_mc->SetMarkerColor(kRed+2);
        hist_mc->SetLineColor(kRed+2);
        hist_signal_subtracted->SetMarkerStyle(20);
        hist_signal_subtracted->SetMarkerSize(0.9);
        hist_signal_subtracted->SetMarkerColor(kBlue+2);
        //---------zoom
        /*hist_mc_zoom->SetMarkerStyle(20);
        hist_mc_zoom->SetMarkerSize(0.9);
        hist_mc_zoom->SetMarkerColor(kRed+2);
        hist_mc_zoom->SetLineColor(kRed+2);
        hist_signal_subtracted_zoom->SetMarkerStyle(20);
        hist_signal_subtracted_zoom->SetMarkerSize(0.9);
        hist_signal_subtracted_zoom->SetMarkerColor(kBlue+2);*/
        //-----------------
        Double_t max_val = TMath::Max(hist_signal_subtracted->GetMaximum(), hist_mc->GetMaximum());
        hist_mc->SetMaximum(max_val * 1.2);
        hist_signal_subtracted->SetMaximum(max_val * 1.2);
        hist_mc->SetMinimum(0);
        hist_signal_subtracted->SetMinimum(0);
        hist_signal_subtracted->Draw("E1X0");        
        hist_mc->Draw("E1X0 SAME");
        hist_signal_subtracted->Draw("E1");
        hist_mc->Draw("E1 SAME");
    
        if(var_mc == "Bchi2cl"){// the label was overlapping with the distribution
            TLegend* legend_compare = new TLegend(0.25, 0.7, 0.68, 0.88);
            legend_compare->SetTextFont(42);
            legend_compare->SetTextSize(0.03);
            legend_compare->SetBorderSize(0);
            legend_compare->SetFillStyle(0);
            legend_compare->AddEntry(hist_mc, "MC", "lep");
            legend_compare->AddEntry(hist_signal_subtracted, "Sideband Subtracted", "lep");
            legend_compare->Draw();}
        else{ 
            TLegend* legend_compare = new TLegend(0.55, 0.7, 0.88, 0.88);
            legend_compare->SetTextFont(42);
            legend_compare->SetTextSize(0.03);
            legend_compare->SetBorderSize(0);
            legend_compare->SetFillStyle(0);
            legend_compare->AddEntry(hist_mc, "MC", "lep");
            legend_compare->AddEntry(hist_signal_subtracted, "Sideband Subtracted", "lep");
            legend_compare->Draw();}

    //if(var_mc=="Bnorm_svpvDistance_2D"){ c_compare->SaveAs(Form("%s_sideband_subtracted_zoom.pdf", var.c_str()));}
    //else{
    c_compare->SaveAs(Form("%s_sideband_subtracted_zoom.pdf", var.c_str()));
    //c_compare->SaveAs(Form("%s_sideband_subtracted.png", var.c_str()));}

/*
    TCanvas* c = new TCanvas(Form("signal_sideband_%s", var_mc.c_str()), var_mc.c_str(), 800, 600); //to check if the subtraction is done properly
    c->SetGrid();
    c->SetTicks();
    gPad->Update();
        hist_signal->SetStats(0);
        hist_sideband->SetStats(0);
        hist_signal_subtracted->SetStats(0);
        hist_signal->SetLineColor(kOrange+7);
        hist_signal->SetFillColor(kOrange+7);
        hist_signal->SetFillStyle(1001);
        hist_sideband->SetLineColor(kBlue+1);
        hist_sideband->SetFillColor(kBlue+1);
        hist_sideband->SetFillStyle(3358);
        hist_signal_subtracted->SetLineColor(kRed+1);
        hist_signal_subtracted->SetFillColorAlpha(kRed+1, 0.35);
        gPad->Update();
        if (hist_signal->Integral() > 0) hist_signal->Scale(1.0 / hist_signal->Integral());
        if (hist_sideband->Integral() > 0) hist_sideband->Scale(1.0 / hist_sideband->Integral());
        hist_signal_subtracted->Scale(1.0 / hist_signal_subtracted->Integral());
        hist_signal->SetTitle(Form(";%s; Entries",var_mc.c_str()));

        Double_t max_val1 = TMath::Max(hist_signal->GetMaximum(), hist_sideband->GetMaximum());
        Double_t max_val2 = TMath::Max(max_val1, hist_signal_subtracted->GetMaximum());

        hist_signal->SetMaximum(max_val2 * 1.2);
        hist_sideband->SetMaximum(max_val2 * 1.2);
        hist_signal->SetMinimum(0);
        hist_sideband->SetMinimum(0);
        hist_signal_subtracted->SetMinimum(0);
        hist_signal->SetTitle("B^{+}");
        hist_signal->Draw("HIST");
        hist_signal_subtracted->Draw("HIST SAME");
        hist_sideband->Draw("HIST SAME");
        
        if(var == "Bchi2cl"){// the label was overlapping with the distribution
            TLegend* legend = new TLegend(0.25, 0.7, 0.68, 0.88);
            legend->SetTextFont(42);
            legend->SetTextSize(0.03);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hist_signal, "Signal Region", "f");
            legend->AddEntry((TObject*)0, Form("signal entries  %0.f", Nsig), "");
            legend->AddEntry(hist_sideband, "Sideband Region", "f");
            legend->AddEntry((TObject*)0, Form("sideband entries %0.f", Nsideband), "");
            legend->AddEntry(hist_signal_subtracted, "After Subtraction","f");
            legend->AddEntry((TObject*)0, Form("subtraction entries %0.f", Nsubtracted), "");
            legend->Draw();}
        
        else{TLegend* legend = new TLegend(0.55, 0.7, 0.88, 0.88);
            legend->SetTextFont(42);
            legend->SetTextSize(0.03);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hist_signal, "Signal Region", "f");
            legend->AddEntry((TObject*)0, Form("signal entries  %0.f", Nsig), "");
            legend->AddEntry(hist_sideband, "Sideband Region", "f");
            legend->AddEntry((TObject*)0, Form("sideband entries %0.f", Nsideband), "");
            legend->AddEntry(hist_signal_subtracted, "After Subtraction","f");
            legend->AddEntry((TObject*)0, Form("subtraction entries %0.f", Nsubtracted), "");
            legend->Draw();}
    c->SaveAs(Form("%s_signal_sideband.pdf", var.c_str()));*/
        //c->SaveAs(Form("%s_signal_sideband.png", var.c_str()));

delete c_compare;
//delete c;
delete hist_mc;
delete hist_signal_subtracted;
delete hist_signal;
delete hist_sideband;

}
}

int main() {
    sbs_var();
}