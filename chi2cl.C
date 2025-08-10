#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <TStyle.h>

#include "ACCSEL.h"

void chi2cl() {



// MC FILEs
const char * files[] = {//MC FILEs
   "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bu_phat5_Bfinder.root",    //ppRef                         
   "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bd_phat5_Bfinder.root",       //ppRef
    "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bs_phat5_Bfinder.root"       //ppRef
   
  // "/lstore/cms/henrique/X3872/MC_DATA/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder.root"
    
           
};
//const char* filesx={  "/lstore/cms/henrique/X3872/MC_DATA/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder.root" };//MC FILEs X3872

const char* files_data={//REAL DATA
     "/lstore/cms/hlegoinha/Bmeson/MC_DATA/DATA_ppref_Bmeson/DATA_ppref_Bmeson.root"  //Bmesons
    //X3872
     //"/lstore/cms/henrique/X3872/MC_DATA/DATA_ppRef_X3872.root"
};

//VARIABLES
int SELplots = 1; //mudar para 1 com ruído e descomentar a linha acima
const char * variables[] = {"Bchi2cl"};
const double ranges[][2] = {{0,1}};
//VARIABLES

/////////////////////////////////  ///////////////////////////  ////////////////

TString cutlevel = ""; // "_RAW", "_ACC", "_SEL", "_TRG", "", 

/////////////////////////////////  ///////////////////////////  ///////////

TString path_to_file = "";

const int nVars = sizeof(variables)/sizeof(variables[0]);

for (int ifile = 0; ifile < sizeof(files)/sizeof(files[0]); ++ifile) {
    path_to_file = Form("%s", files[ifile]);
    //path_to_file = Form("/eos/user/h/hmarques/MC_ppRef_Bmeson/MC_ppRef_Bmeson/%s_Bfinder.root", files[ifile]);

     

    TFile *file = TFile::Open(path_to_file.Data());
    TFile *file_data = TFile::Open(files_data);
    //TFile*filex = TFile::Open(filesx);
    // Get the trees from the file
    TTree *treeMix;
    TTree *treedata;
    //TTree *treex;

    if (path_to_file.Contains("Bs")){                             //Bs
        file->GetObject("Bfinder/ntphi", treeMix);
        file_data->GetObject("Bfinder/ntphi", treedata);
    } else if (path_to_file.Contains("Bd")){                      //Bd
        file->GetObject("Bfinder/ntKstar", treeMix);
        file_data->GetObject("Bfinder/ntKstar", treedata);
    } else if(path_to_file.Contains("Bu")){                       //Bu
        file->GetObject("Bfinder/ntKp", treeMix);
        file_data->GetObject("Bfinder/ntKp", treedata);
    }else{                                                        //X3872
         file->GetObject("Bfinder/ntmix", treeMix);//PSI2S  
        // filex->GetObject("Bfinder/ntmix", treex);//X3872                                       
         file_data->GetObject("Bfinder/ntmix", treedata);
    }


    std::cout << "\n" << "Entries in treeMix: " << treeMix->GetEntries() << std::endl;
    std::cout << "\n" << "Entries in treedata: " << treedata->GetEntries() << std::endl;
    //std::cout << "\n" << "Entries in treex: " << treex->GetEntries() << std::endl;
    for (int i = 0; i < nVars; ++i) {
        TString var = variables[i];

        if(path_to_file.Contains("Bu") && ((var.Contains("trk2") || var.Contains("Ptimb")))) continue; // B+ has less one track!

        // Create a canvas to draw the histograms
        TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetTopMargin(0.05);
        canvas->SetRightMargin(0.05); 

        double hist_Xhigh      = ranges[i][1];
        double hist_Xlow       = ranges[i][0];
        int hist_Nbins;
        int hist_Nbin; 
        //if(var=="Bmass"){hist_Nbin = 150;}  

        ////BINS FOR ZOOM
        if(path_to_file.Contains("Bu")){
            hist_Nbins = 2000;
            hist_Nbin=200;
        }
        else if(path_to_file.Contains("Bd")){
            hist_Nbins = 2000;
            hist_Nbin=200;
        }
        else if(path_to_file.Contains("Bs")){
            hist_Nbins = 2000;
            hist_Nbin=200;
        }
        /////////////////////////
        
        if (var == "nSelectedChargedTracks") {
            hist_Nbin = hist_Xhigh - hist_Xlow;
        } 
        double bin_length_MEV  = (hist_Xhigh - hist_Xlow) / hist_Nbin;
        double bin_length_MEV_Z  = (hist_Xhigh - hist_Xlow) / hist_Nbins;

        if(SELplots){ 
            hist_Nbin=200;  
            if(path_to_file.Contains("Bu")){
                 hist_Nbins = 2000;
            }
            else if(path_to_file.Contains("Bd")){
               hist_Nbins = 2000;
            }
            else if(path_to_file.Contains("Bs")){
                 hist_Nbins = 200;
            }}
        /////////////////////////
       
        
        TString Xlabel ;
        if (var == "Bmass"){ 
            if (path_to_file.Contains("Bs")){
                Xlabel = "m_{J/#Psi K^{+} K^{-}} [GeV/c^{2}]";
            } else if (path_to_file.Contains("Bd")){
                Xlabel = "m_{J/#Psi K^{+} #pi^{-}} [GeV/c^{2}]";
            } else if(path_to_file.Contains("PSI2S")){
                Xlabel="m_{J/#Psi#pi^{+}#pi^{-}} [GeV/c^{2}]";
            } else {
                Xlabel = "m_{J/#Psi K^{+}} [GeV/c^{2}]";
            }
        } else if (var == "Bpt"){ 
            Xlabel = "p_{T} [GeV/c]";
        } else { 
            Xlabel = var.Data();
        }

        // Create histograms
        TH1F *hist_SIG = new TH1F("hist_SIG"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        //TH1F *hist_sig = new TH1F("hist_sig"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh); 
        TH1F *hist_BKG = new TH1F("hist_BKG"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        TH1F *hist     = new TH1F("hist"          , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        TH1F *hist_SIG_WT   = new TH1F("hist_SIG_WT"  , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);        
        TH1F *hist_SIG_BOTH = new TH1F("hist_SIG_BOTH", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        /////////////PLOT ZOOM//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        TH1F *hist_SIGZ= new TH1F("hist_SIGZ"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV_Z) , hist_Nbins, hist_Xlow ,hist_Xhigh);
        TH1F *hist_SIG_WTZ = new TH1F("hist_SIG_WTZ"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV_Z) , hist_Nbins, hist_Xlow ,hist_Xhigh); 
        TH1F *hist_SIG_BOTHZ = new  TH1F("hist_SIG_BOTHZ", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV_Z) , hist_Nbins, hist_Xlow ,hist_Xhigh);
        TH1F *hist_BKGZ = new TH1F("hist_BKGZ"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV_Z) , hist_Nbins, hist_Xlow ,hist_Xhigh);  


        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //SELECT THE acc + presel CUT 

        TString dirNAME = "presel_STUDY/";
        TString Final = "1";      
        TString trgmatches = TRGmatching.Data();   //TRG matching only in ppRef
        TString ACCcuts = "" ;
        TString SELcuts = "" ;

        if (path_to_file.Contains("Bu")){
            ACCcuts    = ACCcuts_ppRef_Bu.Data(); //ppRef
            SELcuts    = SELcuts_ppRef_Bu.Data(); //ppRef
            if (path_to_file.Contains("PbPb")) { 
                ACCcuts = ACCcuts_PbPb_Bu.Data();
                SELcuts = SELcuts_PbPb_Bu.Data();
                trgmatches = "1";
            }
        }
        else {
            ACCcuts    = ACCcuts_ppRef.Data(); //ppRef
            SELcuts    = SELcuts_ppRef.Data(); //ppRef
            if (path_to_file.Contains("PbPb")) { 
                ACCcuts = ACCcuts_PbPb.Data();
                SELcuts = SELcuts_PbPb.Data();
                trgmatches = "1";
            }
        }

        TString cut = "";
        if (cutlevel == "_RAW")       {cut = Form(" %s "                   ,FIDreg.Data());}                                                              //RAW (inside fid reg only)
        else if (cutlevel == "_ACC")  {cut = Form(" %s && %s "             ,FIDreg.Data(), ACCcuts.Data());}                                              //ACC
        else if (cutlevel == "_SEL")  {cut = Form(" %s && %s && %s "       ,FIDreg.Data(), ACCcuts.Data(), SELcuts.Data());}                              //SEL
        else if (cutlevel == "_TRG")  {cut = Form(" %s && %s && %s && %s " ,FIDreg.Data(), ACCcuts.Data(), SELcuts.Data(), trgmatches.Data());}           //TRG
        else if (cutlevel == ""){
            if (!SELplots) {dirNAME  = "";}
            cut = Form(" %s && %s && %s && %s", ACCcuts.Data(), SELcuts.Data(), trgmatches.Data(), Final.Data());                   //Final
        }
        else{
            std::cerr << "Invalid cut level specified: " << cutlevel << std::endl;
            return;
        }                                                                                                 


        TString sepcCASES = "1";
        if (path_to_file.Contains("Bs")){ 
            sepcCASES = "abs(Btktkmass - 1.019455) < 0.015"; // phi meson mass cut 
            treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            treeMix->Draw(Form("%s >> hist_SIGZ", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form(" %s && %s", cut.Data(), sepcCASES.Data()));//((Bmass < 5.289) || (Bmass > 5.449))
            treedata->Draw(Form("%s >> hist_BKGZ", var.Data()), Form(" %s && %s", cut.Data(), sepcCASES.Data()));
        } else if (path_to_file.Contains("Bd")){ 
            sepcCASES = "abs(Btktkmass - 0.89594) < 0.25"; // Kstar meson mass cut
        }

         
         //treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));  // SIG
        if (path_to_file.Contains("Bd")){ 
            if(var== "Bmass"){treedata->Draw(Form("%s >> hist_BKG"     , var.Data()), Form(" %s && %s", cut.Data(), sepcCASES.Data()));}//  && ((Bnorm_svpvDistance_2D > 16) && (Bnorm_svpvDistance_2D < 72)) && (Btrk2dR<0.6)  &&(Bnorm_svpvDistance_2D>3.9916)  && 
            else{ // (Btrk2dR<0.6)  &&
                treeMix->Draw(Form("%s >> hist_SIG_WT"  , var.Data()), Form(" (Bgen == 41000)  && Bnorm_svpvDistance_2D>2 && %s && %s", cut.Data(), sepcCASES.Data()));            
                treeMix->Draw(Form("%s >> hist_SIG_WTZ"  , var.Data()), Form(" (Bgen == 41000)  && Bnorm_svpvDistance_2D>2 && %s && %s", cut.Data(), sepcCASES.Data()));                    // WT component && (Balpha<0.01) && ((Bnorm_svpvDistance_2D > 16) && (Bnorm_svpvDistance_2D < 72)) 
	            treedata->Draw(Form("%s >> hist_BKG"     , var.Data()), Form(" ((Bmass < 5.17807 ) || (Bmass >5.38111))&&   %s && %s", cut.Data(), sepcCASES.Data()));// (Balpha<0.01)  &&  ((Bnorm_svpvDistance_2D > 16) && (Bnorm_svpvDistance_2D < 72))&&
                treedata->Draw(Form("%s >> hist_BKGZ"     , var.Data()), Form(" ((Bmass < 5.17807 ) || (Bmass >5.38111))&&   %s && %s", cut.Data(), sepcCASES.Data()));
                treeMix->Draw(Form("%s >> hist_SIG_BOTH", var.Data()), Form("(%s || (Bgen == 41000))  && Bnorm_svpvDistance_2D>2 && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));  // SIG + WT  && (Balpha<0.01) &&  ((Bnorm_svpvDistance_2D > 16) && (Bnorm_svpvDistance_2D < 72)) && Bnorm_svpvDistance_2D>3.9916
                treeMix->Draw(Form("%s >> hist_SIG_BOTHZ", var.Data()), Form("(%s || (Bgen == 41000))  && Bnorm_svpvDistance_2D>2 && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data())); 
            }
        } else if(path_to_file.Contains("Bu")) {
            if(var=="Bmass"){
                treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form(" (Balpha<0.015) && %s && %s", cut.Data(), sepcCASES.Data()));
            } 
            else{
                treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form(" Bnorm_svpvDistance_2D>2 &&%s && %s && %s ", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));  // SIG  (Btrk1dR<1.285)  && (Balpha<0.015) 
                treeMix->Draw(Form("%s >> hist_SIGZ", var.Data()), Form(" Bnorm_svpvDistance_2D>2 &&%s && %s && %s ", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
                treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form(" Bmass>5.38772 &&  %s && %s", cut.Data(), sepcCASES.Data()));  // (Bmass>5.17501) && (Balpha<0.015)
                treedata->Draw(Form("%s >> hist_BKGZ", var.Data()), Form(" Bmass>5.38772 &&  %s && %s", cut.Data(), sepcCASES.Data())); 
            }}
 
        //SELECT THE acc + presel CUT 
        
        // Customize the Histograms
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);

        hist_SIG->SetLineColor(kOrange+7);
        hist_SIG->SetFillColor(kOrange+7); 
        hist_SIGZ->SetLineColor(kOrange+7);
        hist_SIGZ->SetFillColor(kOrange+7); 

        //hist_SIG->SetFillStyle(3001); 
        //hist_SIG->SetLineWidth(2);
        //hist_SIG->SetLineStyle(2);

        hist_SIG_BOTH->SetLineColor(kOrange+7);
        hist_SIG_BOTH->SetFillColor(kOrange+7);
        hist_SIG_BOTHZ->SetLineColor(kOrange+7);
        hist_SIG_BOTHZ->SetFillColor(kOrange+7);
        hist_SIG_WT->SetLineColor(kOrange);
        hist_SIG_WT->SetFillColor(kOrange);
        hist_SIG_WTZ->SetLineColor(kOrange);
        hist_SIG_WTZ->SetFillColor(kOrange);    

        hist_BKG->SetLineColor(kBlue);
        hist_BKG->SetFillColor(kBlue);     
        hist_BKG->SetFillStyle(3358);
        hist_BKG->SetMaximum(1.1 * hist_BKG->GetMaximum()); 
        hist_BKGZ->SetLineColor(kBlue);
        hist_BKGZ->SetFillColor(kBlue);     
        hist_BKGZ->SetFillStyle(3358);
        //hist_BKGZ->SetMaximum(1.1 * hist_BKGZ->GetMaximum());// Increase the max range to give some space
        //hist_BKG->SetLineStyle(2);
        //hist_BKG->SetLineWidth(2);

        //hist_sig->SetLineColor(kOrange+7);
        //hist_sig->SetFillColor(kOrange+7);
        //hist_sig->SetFillStyle(3001);
        //hist_BKG->SetMarkerStyle(20); // Circle marker
        //hist_BKG->SetMarkerSize(.8); // Bigger dots
        // Customize the Histograms

        hist_SIG->SetName("MC_SIG");  // <--- This affects the stat box label
        hist_SIGZ->SetName("MC_SIG");
        hist_BKG->SetName("DATA_BKG");  // <--- Also affects the stat box
        hist_BKGZ->SetName("DATA_BKG");
        hist_SIG_BOTH->SetName("MC_SIG");
        hist_SIG_BOTHZ->SetName("MC_SIG");
       /* if (path_to_file.Contains("Bu")) {
            hist_BKG->SetName("DATA ");
        } */


        if(path_to_file.Contains("X3872")){
            //hist_sig->SetName("MC_SIG_X3872");
            hist_SIG->SetName("MC_SIG_PSI2S");
                }
        if (SELplots == 1) { // NORMALIZE
              double int_sig     = hist_SIG->Integral();
              double int_sigZ   = hist_SIGZ->Integral();
             // double int_sig_x   = hist_sig->Integral();
              double int_bkg     = hist_BKG->Integral();
              double int_bkgZ   = hist_BKGZ->Integral();
              double int_sig_wt  = hist_SIG_BOTH->Integral();
              double int_sig_wtZ = hist_SIG_BOTHZ->Integral();

            if (int_sig > 0|| int_sig_wt > 0 ||int_bkg > 0){
                hist_SIG->Scale(1.0 / int_sig);
                hist_SIGZ->Scale(1.0 / int_sigZ);
                hist_BKG->Scale(1.0 / int_bkg);
                hist_BKGZ->Scale(1.0 / int_bkgZ);
                hist_SIG_BOTH->Scale(1.0 / int_sig_wt);
                hist_SIG_BOTHZ->Scale(1.0 / int_sig_wtZ);
               // hist_SIG_WT->Scale(1.0 / int_sig_wt);
               // hist_SIG_WTZ->Scale(1.0 / int_sig_wtZ);
            
          }
        } 

        if(1){// set the y-axis maximum if needed
            Double_t     max_val = TMath::Max(hist->GetMaximum(), TMath::Max(hist_BKG->GetMaximum(), hist_SIG->GetMaximum()));
            Double_t max_val_f = TMath::Max(hist->GetMaximum(), TMath::Max(hist_BKGZ->GetMaximum(), hist_SIGZ->GetMaximum()));
            if(SELplots) {
                if (path_to_file.Contains("Bd")) {
                    max_val =TMath::Max( hist_BKG->GetMaximum(), hist_SIG_BOTH->GetMaximum()) ;
                    max_val_f = TMath::Max(hist_BKGZ->GetMaximum(), hist_SIG_BOTHZ->GetMaximum());
                    hist_SIG_BOTH->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
                    hist_SIG_BOTHZ->SetMaximum(hist_SIG_BOTHZ->GetMaximum() * 1.1);  // Increase the max range to give some space
                    //hist_BKG->SetMaximum(max_val * 1.1);
                } else {
                    max_val = TMath::Max( hist_SIG->GetMaximum(), hist_BKG->GetMaximum());
                    max_val_f = TMath::Max(hist_SIGZ->GetMaximum(), hist_BKGZ->GetMaximum());
                    hist_SIG->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
                    hist_SIGZ->SetMaximum(hist_SIGZ->GetMaximum()* 1.1);  // Increase the max range to give some space
                    hist_BKG->SetMaximum(max_val * 1.1); 
                    hist_BKGZ->SetMaximum(hist_SIGZ->GetMaximum()* 1.1);  // Increase the max range to give some space
            } 
        }}

        // Draw the histograms
            hist->SetStats(0);

        if (SELplots && path_to_file.Contains("Bd")){
            if (var == "Bmass"){
                hist_BKG->Draw("HIST ");
            }else{
                hist_SIG_BOTH->Draw("HIST");
                //hist_SIG_WT->Draw("HIST");
                hist_BKG->Draw("HIST SAMES");
            }} 
        else if(SELplots && path_to_file.Contains("Bu")){
             if (var == "Bmass"){
                hist_BKG->Draw("HIST ");
            }else{
                hist_SIG->Draw("HIST");
                hist_BKG->Draw("HIST SAMES");
            }} 
        else if(SELplots && path_to_file.Contains("Bs")){
            if(var=="Bmass"){
                hist_BKG->Draw("HIST");
            }else{
                hist_BKG->Draw("HIST");
                hist_SIG->Draw("HIST SAMES");}}

               
                
        if(!SELplots) hist->Draw("HIST SAME");
        canvas->Update();
        gPad->Update();


        TPad *subpad = new TPad("subpad","subpad",0.4,0.15,0.9,0.65);
       /* subpad->SetFrameLineWidth(2);
        subpad->SetFrameLineColor(kBlack);
        subpad->SetFrameLineStyle(1);  */ // solid line
        // after creating / cd’ing into the sub-pad
          // ticks on right edge
        subpad->SetTopMargin(0.05);
        subpad->SetBottomMargin(0.2);
        subpad->SetLeftMargin(0.2);
        subpad->SetRightMargin(0.05);
        subpad->Draw();
        subpad->cd();
        
        // Disable the statistics box
        if(path_to_file.Contains("Bu")){
            gPad->SetTickx(1);   // ticks on top edge
            gPad->SetTicky(1); 
            TH1D *h2=(TH1D *)hist_SIGZ->Clone("hist");
            TH1D *h2_bkg=(TH1D *)hist_BKGZ->Clone("hist same");
            h2->SetStats(0);
            h2_bkg->SetStats(0);
            h2->Draw("hist");
            h2_bkg->Draw("hist same");
            TAxis *xaxis=h2->GetXaxis();
            TAxis *yaxis=h2->GetYaxis();
            xaxis->SetRangeUser(0,0.05);
            xaxis->SetTitle(Xlabel.Data());
            yaxis->SetTitle("Normalized Entries"); 
            h2->GetXaxis()->SetTitleSize(0.06);
            h2->GetYaxis()->SetTitleSize(0.06);   
            h2->GetXaxis()->SetLabelSize(0.04);   // 5 % of pad height
            h2->GetYaxis()->SetLabelSize(0.05);   // 5 % of pad width
            xaxis->SetTickSize(0.04);
            yaxis->SetTickSize(0.04);}
        else if(path_to_file.Contains("Bd")){
            gPad->SetTickx(1);   // ticks on top edge
            gPad->SetTicky(1); 
            TH1D *h2=(TH1D *)hist_SIG_BOTHZ->Clone("hist");
            //TH1D *h2_sigWT=(TH1D *)hist_SIG_WT->DrawCopy("hist same");
            TH1D *h2_bkg=(TH1D *)hist_BKGZ->Clone("hist same");
            h2->SetStats(0);
            h2_bkg->SetStats(0);
            h2->Draw("hist");
            h2_bkg->Draw("hist same");
            TAxis *xaxis=h2->GetXaxis();
            TAxis *yaxis=h2->GetYaxis();
            xaxis->SetRangeUser(0,0.05);
            //yaxis->SetRangeUser(0,0.02);
            xaxis->SetTitle(Xlabel.Data());
            yaxis->SetTitle("Normalized Entries"); 
            h2->GetXaxis()->SetTitleSize(0.06);
            h2->GetYaxis()->SetTitleSize(0.06);   
            h2->GetXaxis()->SetLabelSize(0.04);   // 5 % of pad height
            h2->GetYaxis()->SetLabelSize(0.05);   // 5 % of pad width
            xaxis->SetTickSize(0.04);
            yaxis->SetTickSize(0.04);}
            else if(path_to_file.Contains("Bs")){
            gPad->SetTickx(1);   // ticks on top edge
            gPad->SetTicky(1); 
            TH1D *h2=(TH1D *)hist_SIGZ->DrawCopy("hist");
            TH1D *h2_bkg=(TH1D *)hist_BKGZ->DrawCopy("hist same");
            TAxis *xaxis=h2->GetXaxis();
            TAxis *yaxis=h2->GetYaxis();
            xaxis->SetRangeUser(0,0.05);
            //yaxis->SetRangeUser(0,0.025);
            xaxis->SetTitle(Xlabel.Data());
            yaxis->SetTitle("Normalized Entries"); 
            h2->GetXaxis()->SetTitleSize(0.06);
            h2->GetYaxis()->SetTitleSize(0.06);   
            h2->GetXaxis()->SetLabelSize(0.05);   // 5 % of pad height
            h2->GetYaxis()->SetLabelSize(0.05);   // 5 % of pad width
            xaxis->SetTickSize(0.04);
            yaxis->SetTickSize(0.04);}
    

        canvas->cd();
        canvas->Update();
        gPad->Update();

        // Move and color the stat boxes
        TPaveStats *st_bkg = (TPaveStats*)hist_BKG->GetListOfFunctions()->FindObject("stats");
        if (st_bkg) {
            st_bkg->SetTextColor(kBlue);
            st_bkg->SetLineColor(kBlue); 
            st_bkg->SetX1NDC(0.75);
            st_bkg->SetX2NDC(0.95);
            st_bkg->SetY1NDC(0.85);
            st_bkg->SetY2NDC(0.95);
            st_bkg->Draw();
        }

        TPaveStats *st_bkgz = (TPaveStats*)hist_BKGZ->GetListOfFunctions()->FindObject("stats");
        if (st_bkgz) {
            st_bkgz->SetTextColor(kBlue);
            st_bkgz->SetLineColor(kBlue); 
            st_bkgz->SetX1NDC(0.75);
            st_bkgz->SetX2NDC(0.95);
            st_bkgz->SetY1NDC(0.85);
            st_bkgz->SetY2NDC(0.95);
            st_bkgz->Draw();
        }

       TPaveStats *st_sig = (TPaveStats*)hist_SIG->GetListOfFunctions()->FindObject("stats");
        if (st_sig) {
            st_sig->SetTextColor(kOrange+7);
            st_sig->SetLineColor(kOrange+7);
            st_sig->SetX1NDC(0.75);
            st_sig->SetX2NDC(0.95);
            st_sig->SetY1NDC(0.75);
            st_sig->SetY2NDC(0.85);
            st_sig->Draw();
        }

        TPaveStats *st_sigz = (TPaveStats*)hist_SIGZ->GetListOfFunctions()->FindObject("stats");
        if (st_sigz) {
            //st_sig->SetTextColor(kOrange+7);
            //st_sig->SetLineColor(kOrange+7);
            //st_sig->SetLineColor(kRed);         visible 1
            st_sigz->SetTextColor((kRed+1));
            st_sigz->SetLineColor((kRed+1));
            st_sigz->SetX1NDC(0.75);             
            st_sigz->SetX2NDC(0.95);
            st_sigz->SetY1NDC(0.75);
            st_sigz->SetY2NDC(0.85);
            st_sigz->Draw();
        }
    
        TPaveStats *st_sigWT = (TPaveStats*)hist_SIG_WT->GetListOfFunctions()->FindObject("stats");
        if (st_sigWT) {
            st_sigWT->SetTextColor(kOrange);
            st_sigWT->SetLineColor(kOrange);
            st_sigWT->SetX1NDC(0.75);
            st_sigWT->SetX2NDC(0.95);
            st_sigWT->SetY1NDC(0.65);
            st_sigWT->SetY2NDC(0.75);
            st_sigWT->Draw();
        }     
        TPaveStats *st_sigWTZ = (TPaveStats*)hist_SIG_WTZ->GetListOfFunctions()->FindObject("stats");
        if (st_sigWTZ) {
            st_sigWTZ->SetTextColor(kOrange);
            st_sigWTZ->SetLineColor(kOrange);
            st_sigWTZ->SetX1NDC(0.75);
            st_sigWTZ->SetX2NDC(0.95);
            st_sigWTZ->SetY1NDC(0.65);
            st_sigWTZ->SetY2NDC(0.75);
            st_sigWTZ->Draw();
        }    
         TPaveStats *st_sig_both = (TPaveStats*)hist_SIG_BOTH->GetListOfFunctions()->FindObject("stats");
        if (st_sig_both) {
            st_sig_both->SetTextColor(kOrange+7);
            st_sig_both->SetLineColor(kOrange+7);
            st_sig_both->SetX1NDC(0.75);
            st_sig_both->SetX2NDC(0.95);
            st_sig_both->SetY1NDC(0.75);
            st_sig_both->SetY2NDC(0.85);
            st_sig_both->Draw();
        }    

        TPaveStats *st_sig_bothz = (TPaveStats*)hist_SIG_BOTHZ->GetListOfFunctions()->FindObject("stats");
        if (st_sig_bothz) {
            st_sig_bothz->SetTextColor(kOrange+7);
            st_sig_bothz->SetLineColor(kOrange+7);
            st_sig_bothz->SetX1NDC(0.75);
            st_sig_bothz->SetX2NDC(0.95);
            st_sig_bothz->SetY1NDC(0.75);
            st_sig_bothz->SetY2NDC(0.85);
            st_sig_bothz->Draw();
        } 
    

        // LATEX text
        if(0){
            double Nsignal = hist_SIG->GetEntries();
            double Nbkg = hist_BKG->GetEntries();
            double significance = (Nbkg > 0) ? Nsignal / sqrt(Nbkg) : 0;

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.022);
            latex.SetTextColor(kOrange+7); // Same as hist_SIG
            latex.DrawLatex(0.18, 0.82, Form("N_{sig} = %.0f", Nsignal));
            latex.SetTextColor(kBlue);     // Same as hist_BKG
            latex.DrawLatex(0.18, 0.85, Form("N_{bkg} = %.0f", Nbkg));
        }
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Add a legend
        auto legend = new TLegend(0.15, 0.7, 0.25, 0.9);
        legend->AddEntry(hist_SIG, "MC SIG", "l");
        legend->AddEntry(hist_BKG, "MC BKG", "l");

        //legend->Draw();
              


        // Save the canvas as an image
        TString particleNAME = "Bu";
        TString systemNAME = "MC_ppRef_";
        if (path_to_file.Contains("Bs")){
            particleNAME = "Bs";
        } else if (path_to_file.Contains("Bd")){
            particleNAME = "Bd";
        } 
        else if(path_to_file.Contains("Rho")){
            particleNAME="X3872";
        }else if(path_to_file.Contains("PSI2S")){
            particleNAME="PSI2S";}
        if (path_to_file.Contains("PbPb"))  { systemNAME = "MC_PbPb_";}

        canvas->SaveAs(Form("./%s%s%s_%s%s_2000bins_200bins.pdf", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));
        //canvas->SaveAs(Form("./%s%s%s_%s%s_optimalcut.root", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));


        // Clean up
        delete hist_SIG;
        delete hist_SIGZ;
        delete hist_SIG_WT;
        delete hist_SIG_WTZ;
        delete hist_SIG_BOTH;
        delete hist_SIG_BOTHZ;
        delete hist_BKG;
        delete hist_BKGZ;
        delete hist;
        delete canvas;
        
    }

}
}
    

int main() {
    chi2cl();
    return 0;
}