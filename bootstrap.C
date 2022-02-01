void bootstrap(){
  const Int_t nOfSamples = 10;

  TFile* fileM = TFile::Open("plotsWithEfficiency/pPb/assRef/v2_TF_Proton_sample_0.root");
  fileM->ls();
  // TH1D* histM = (TH1D*)fileM->Get("hFinalV2_Proton_6");
  TH1D* histM = (TH1D*)fileM->Get("hFinalV2_Proton_1");

  TList* list = new TList();

  TFile* file[nOfSamples];
  TH1D* hist[nOfSamples];
  for(Int_t iSam(0); iSam < nOfSamples; iSam++){
    file[iSam] = TFile::Open(Form("plotsWithEfficiency/pPb/samples/v2_TF_Proton_sample_%d.root",iSam));
    // hist[iSam] = (TH1D*)file[iSam]->Get("hFinalV2_Proton_6");
    hist[iSam] = (TH1D*)file[iSam]->Get("hFinalV2_Proton_1");
    list->Add(hist[iSam]);
  }

  TH1D* bootstrap = (TH1D*) histM->Clone("BS");
  bootstrap->Sumw2();






  TRandom *rand = new TRandom();
  std::vector<Double_t> bootstrapVector;
  for(Int_t iBin(0); iBin < histM->GetNbinsX()+2; iBin++){
    Int_t iCount = 0;
    Double_t dValue = 0.0;
    Double_t centralValue = histM->GetBinContent(iBin);

    Double_t dMeanBootStrap = 0.0;
    Double_t dVarBootStrap = 0.0;
    Double_t dCount = 0.0;

    for(Int_t bootSample(0); bootSample < 100; bootSample++){
      Double_t dSum = 0.0;
      Double_t dSumWeights = 0.0;

      //creating random dataset
      for(Int_t combi(0); combi < nOfSamples; combi++){
        Int_t index = (Int_t) rand->Uniform(0,nOfSamples);

        Double_t dContent = hist[index]->GetBinContent(iBin);
        Double_t dError = hist[index]->GetBinError(iBin);
        // Double_t dError = 1.0;
        Double_t dWeight = 0;
        if(dError > 0) dWeight = TMath::Power(1./dError,2.0);
        if(dWeight > 0) {
          dSum += dContent*dWeight;
          dSumWeights += dWeight;
        }
      } // end combi

      if(dSumWeights > 0) {
        Double_t help = dSum/dSumWeights;
        bootstrapVector.push_back(help);
        dMeanBootStrap += help;
        dCount++;
      }
    } //end bootSample

    if(dCount > 0) dMeanBootStrap = dMeanBootStrap/dCount;
    Int_t bootSize = bootstrapVector.size();
    for(Int_t bootSample(0); bootSample < bootSize; bootSample++){
      Double_t help = bootstrapVector[bootSample] - dMeanBootStrap;
      dVarBootStrap += TMath::Power(help,2);
    } //end bootSample

    dVarBootStrap = dVarBootStrap/(bootSize-1);
    bootstrap->SetBinError(iBin, TMath::Sqrt(dVarBootStrap));
    bootstrapVector.clear();
  } // end bins


  new TCanvas();


  histM->SetLineColor(kRed);
  histM->Draw();
  bootstrap->Draw("same");

  TH1D* ratioUnc = (TH1D*) histM->Clone("ratioUnc");
  for(Int_t iBin(0); iBin < ratioUnc->GetNbinsX()+2; iBin++){
    Double_t uncBS = bootstrap->GetBinError(iBin);
    Double_t uncSta = histM->GetBinError(iBin);
    if(uncBS > 0.0) ratioUnc->SetBinContent(iBin, uncSta / uncBS);
    ratioUnc->SetBinError(iBin,0.0);
  }

  new TCanvas();
  ratioUnc->Draw();



  TFile* fileOut = TFile::Open("BS_TF_Proton.root", "RECREATE");
  ratioUnc->Write();
  bootstrap->Write();
  histM->Write();




  return;
}
