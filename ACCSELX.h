    TString ACCcuts_ppRef = 
        // Trk Acceptance
        "(abs(Btrk1Eta) < 2.4 && abs(Btrk2Eta) < 2.4)"       
        "&& (Btrk1Pt  > 0.5 && Btrk2Pt  > 0.5)"                     // <--------
        // muon Acceptance   
        "&& (Bmu1isAcc == 1 && Bmu2isAcc == 1)";
       // "&& Bchi2cl > 0.05"; // B+ chi2cl cut
    
  

    TString ACCcuts_PbPb = 
        // Trk Acceptance
        "(abs(Btrk1Eta) < 2.4 && abs(Btrk2Eta) < 2.4)"       
        "&& (Btrk1Pt  > 0.9 && Btrk2Pt  > 0.9)"                     // <--------
        // muon Acceptance   
        "&& (Bmu1isAcc == 1 && Bmu2isAcc == 1)"; 
        //"&& Bchi2cl > 0.05"; // B+ chi2cl cut  
    

    TString FIDreg =
        "(Bpt >= 2 )"                    //FID reg pT
        "&& (abs(By)  <= 2.4)";               //FID reg y

    TString SELcuts_ppRef =
        //Trk Selection
        "((Btrk1PtErr / Btrk1Pt) < 0.1 && (Btrk2PtErr / Btrk2Pt) < 0.1)"
        "&& ((Btrk1nPixelLayer + Btrk1nStripLayer) > 10 && (Btrk2nPixelLayer + Btrk2nStripLayer) > 10)"
        "&& ((Btrk1Chi2ndf/(Btrk1nPixelLayer + Btrk1nStripLayer)) < 0.18 && (Btrk2Chi2ndf/(Btrk2nPixelLayer + Btrk2nStripLayer)) < 0.18)"
        "&& Btrk2highPurity && Btrk1highPurity"
        // Muon Selection
        "&& (Bmu1SoftMuID == 1 && Bmu2SoftMuID == 1)"               // <--------
        // Dimuon Selection
        "&& abs(Bujmass - 3.096916) < 0.15 "
        "&& BujvProb > 0.01";

    TString SELcuts_PbPb =
        //Trk Selection
        "((Btrk1PtErr / Btrk1Pt) < 0.1 && (Btrk2PtErr / Btrk2Pt) < 0.1)"
        "&& ((Btrk1nPixelLayer + Btrk1nStripLayer) > 10 && (Btrk2nPixelLayer + Btrk2nStripLayer) > 10)"
        "&& ((Btrk1Chi2ndf/(Btrk1nPixelLayer + Btrk1nStripLayer)) < 0.18 && (Btrk2Chi2ndf/(Btrk2nPixelLayer + Btrk2nStripLayer)) < 0.18)"
        "&& Btrk2highPurity && Btrk1highPurity"
        // Muon Selection
        "&& (Bmu1HybridSoftMuID == 1 && Bmu1HybridSoftMuID == 1)"   // <--------
        // Dimuon Selection
        "&& abs(Bujmass - 3.096916) < 0.15 "
        "&& BujvProb > 0.01";



        

    TString isMCsignal = "((Bgen == 23333) || (Bgen == 24333) || (Bgen == 23433) || (Bgen == 24433))";

    TString TRGmatching = "(Bmu1isTriggered == 1 && Bmu2isTriggered == 1 )"; // B+ chi2cl cut








////// B+ has less one track!

    TString ACCcuts_ppRef_Bu = 
        // Trk Acceptance
        "(abs(Btrk1Eta) < 2.4) && (Btrk1Pt  > 0.5)"                     // <--------
        // muon Acceptance   
        "&& (Bmu1isAcc == 1 && Bmu2isAcc == 1)" 
        //B mesons selection
        //"&& Bnorm_svpvDistance >2";
        //"&& Bchi2cl > 0.05"; // B+ chi2cl cut

        

    TString ACCcuts_PbPb_Bu = 
        // Trk Acceptance
        "(abs(Btrk1Eta) < 2.4 && (Btrk1Pt  > 0.9)"                     // <--------
        // muon Acceptance   
        "&& (Bmu1isAcc == 1 && Bmu2isAcc == 1)"
        //B mesons selection   
        //"&& Bnorm_svpvDistance >2";
        //"&& Bchi2cl > 0.05"; // B+ chi2cl cut

    TString SELcuts_ppRef_Bu =
        //Trk Selection
        "(Btrk1PtErr / Btrk1Pt) < 0.1"
        "&& (Btrk1nPixelLayer + Btrk1nStripLayer) > 10 "
        "&& (Btrk1Chi2ndf/(Btrk1nPixelLayer + Btrk1nStripLayer)) < 0.18 "
        "&& Btrk1highPurity"
        // Muon Selection
        "&& (Bmu1SoftMuID == 1 && Bmu2SoftMuID == 1)"               // <--------
        // Dimuon Selection
        "&& abs(Bujmass - 3.096916) < 0.15 "
        "&& BujvProb > 0.01";

    TString SELcuts_PbPb_Bu =
        //Trk Selection
        "(Btrk1PtErr / Btrk1Pt) < 0.1"
        "&& (Btrk1nPixelLayer + Btrk1nStripLayer) > 10 "
        "&& (Btrk1Chi2ndf/(Btrk1nPixelLayer + Btrk1nStripLayer)) < 0.18 "
        "&& Btrk1highPurity"
        // Muon Selection
        "&& (Bmu1HybridSoftMuID == 1 && Bmu1HybridSoftMuID == 1)"   // <--------
        // Dimuon Selection
        "&& abs(Bujmass - 3.096916) < 0.15 "
        "&& BujvProb > 0.01";

        
   
  
