#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binst[]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
Double_t binst[]={0.5,0.7,0.9,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.5,0.7,0.9};
// Double_t binst[]={0.5,0.7,0.9,1.25};
// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0};
// Double_t binst[]={0.2,0.5,1.0,1.5,2.0,3.0,4.0};
// Double_t binst[]={0.5,2.0};
// Double_t binst[]={1.8,4.8}; // FMDA
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

const Double_t binsa[]={-3.2,-1.8}; // FMDC
// const Double_t binsa[]={1.8,4.8}; // FMDA
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

Double_t binscmin[] = {0.,60.};
Double_t binscmax[] = {20.,100.};
const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

Bool_t withAss = kFALSE;
Bool_t sameBin = kFALSE;
Bool_t useEff = kFALSE;
Bool_t isPP = kFALSE;

enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSample,kMass,kPt_TPC_trig};
enum {triMinv, triPt, triPVZ, triSample};

void readFMD_V0(TString corr = "TPCFMDC", TString species = "Lambda", Int_t ic = 0, TString cent = "peri", TString mcType = "gen")
{
 Double_t cen1;
 Double_t cen2;
 if(ic == 0)
 {
  cen1 = binscmin[0];
  cen2 = binscmax[0];;
 }

 if(ic == 1)
 {
  cen1 = binscmin[1];
  cen2 = binscmax[1];;
 }


 Double_t binscmin[] = {cen1};
 Double_t binscmax[] = {cen2};
 const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

 const char* strPtTrg;
 if(corr == "FMDAFMDC"){
   strPtTrg = "#eta_{FMD-A}";
 }
 else{
   strPtTrg = "p_{T}^{Tri} (GeV/c)";
 }
 const char* strPtAss;
 if ( corr == "TPCFMDA") {
  strPtAss = "#eta_{FMD-A}";
 }
 else{
   strPtAss = "#eta_{FMD-C}";
 }

 TFile *f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_LHC16q/train_4197/AR_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data()),"READ");
 // TFile *f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_MC/train_1328/AR_%s_%s_%s.root",mcType.Data(),cent.Data(),species.Data()),"READ");
 if(!f) { printf("f\n"); return; }
 AliCFGridSparse *sparSig = nullptr;
 if(species == "Charged") sparSig = (AliCFGridSparse*)f->Get("fhChargedSE_SelStep0");
 else sparSig = (AliCFGridSparse*)f->Get(Form("fhPidSE_%s_SelStep0",species.Data()));
 if(!sparSig) { printf("sparSig\n"); return; }

 // AliCFGridSparse *sparTriSig = (AliCFGridSparse*)f->Get("sparTriSig");
 AliCFGridSparse *sparMix = nullptr;
 if(species == "Charged") sparMix = (AliCFGridSparse*)f->Get("fhChargedME_SelStep0");
 else sparMix = (AliCFGridSparse*)f->Get(Form("fhPidME_%s_SelStep0",species.Data()));
 if(!sparMix) { printf("sparMix\n"); return; }

 AliCFGridSparse *sparTrigger = nullptr;
 if(species == "Charged") sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_SelStep0"));
 else sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));
 if(!sparTrigger) printf("Problem with sparTrigger! \n ");

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = "";
 // outFileName = Form("FMD_Data/Result_%s_%s_%s_test.root",corr.Data(),cent.Data(),species.Data());
 outFileName = Form("FMD_Data/Purity_%s.root",species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");

 Double_t binsV2[14] = {0.5,0.7,0.9,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
 TH1D* purity = new TH1D(Form("purity_%s",species.Data()), "Purity; #it{p}_{T}; Purity", 13, binsV2);

 Int_t nTrig[3] = {0,0,0};

    for (Int_t it=0;it<nt;it++) {
      TH2D* hPhiEtaSMsum[3]={nullptr};

      Double_t minfit = 0.0;
      Double_t maxfit = 0.0;
      if(species == "K0s"){
        minfit = 0.475;
        maxfit = 0.525;
      }
      else{
        minfit = 1.08;
        maxfit = 1.15;
      }

      //first pT range
      sparTrigger->GetAxis(triMinv)->SetRangeUser(minfit,maxfit);
      sparTrigger->GetAxis(triPt)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      TH1D* minv = (TH1D*)sparTrigger->Project(triMinv);
      minv->SetTitle(Form("Projection, M_{inv}, %.2f < p_{T} < %.2f",binst[it],binst[it+1]));

      TF1* fitMinv = new TF1(Form("fitMinv_%d",it), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*([5]*TMath::Gaus(x,[6],[7])+(1.0-[5])*TMath::Gaus(x,[8],[9]))",minfit,maxfit);

      Double_t parDef = 0.0;
      Double_t parLow = 0.0;
      Double_t parHigh = 0.0;
      //mean
      if(species == "K0s"){
        parDef = 0.4976;
        parLow = 0.49;
        parHigh = 0.505;
      }
      else{
        parDef = 1.115;
        parLow = 1.1;
        parHigh = 1.121;
      }
      fitMinv->SetParameter(6, parDef);
      fitMinv->SetParLimits(6, parLow, parHigh);
      fitMinv->SetParameter(8, parDef);
      fitMinv->SetParLimits(8, parLow, parHigh);
      //sigma
      if(species == "K0s"){
        parDef = 0.001;
        parLow = 0.001;
        parHigh = 0.005;
      }
      else{
        parDef = 0.001;
        parLow = 0.001;
        parHigh = 0.02;
      }
      fitMinv->SetParameter(7, parDef);
      fitMinv->SetParLimits(7, parLow, parHigh);
      fitMinv->SetParameter(9, parDef);
      fitMinv->SetParLimits(9, parLow, parHigh);

      minv->Fit(fitMinv, "QMI");

      minv->SetTitle(Form("minv_%d",it));
      minv->Write();
      fitMinv->Write();

      TF1* bkgMinv = new TF1(Form("bkgMinv_%d",it), "[0] + [1]*x + [2]*x*x + [3]*x*x*x",minfit,maxfit);
      bkgMinv->SetParameter(0, fitMinv->GetParameter(0));
      bkgMinv->SetParameter(1, fitMinv->GetParameter(1));
      bkgMinv->SetParameter(2, fitMinv->GetParameter(2));
      bkgMinv->SetParameter(3, fitMinv->GetParameter(3));
      bkgMinv->Write();

      Double_t meanPDG = 0.0;
      if(species == "K0s") meanPDG = 0.497614;
      else meanPDG = 1.11568;

      Double_t sigmaFit = fitMinv->GetParameter(7);
      if(fitMinv->GetParameter(9) > sigmaFit) sigmaFit = fitMinv->GetParameter(9);

      Double_t minRange = meanPDG - 3.0*sigmaFit;
      Double_t maxRange = meanPDG + 3.0*sigmaFit;

      printf("\n range: %f -- %f \n\n", minRange, maxRange);
      printf("\n\n\n sigma1: \t %f \t\t sigma 2: \t %f \n", fitMinv->GetParameter(7), fitMinv->GetParameter(9));

      Double_t signal = fitMinv->Integral(minRange+0.005,maxRange-0.005);
      Double_t background = bkgMinv->Integral(minRange+0.005,maxRange-0.005);
      Double_t pur = (signal)/(signal+background);
      purity->SetBinContent(it+1,pur);
      continue;

      if(withAss && it == 1) continue;

      //pt or dphi trigger loop
      sparSig->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparMix->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      Double_t tmin = sparSig->GetAxis(kPt_TPC_trig)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_trig)->GetFirst());
      Double_t tmax = sparSig->GetAxis(kPt_TPC_trig)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_trig)->GetLast());

      sparTrigger->GetAxis(triPt)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);

      for(Int_t iM(0); iM < 2; iM++){

        if(iM == 1) continue;

        TH1D* hTriggersS = nullptr;
        if(iM == 0){
          sparTrigger->GetAxis(triMinv)->SetRangeUser(minRange+0.005,maxRange-0.005);
          hTriggersS = (TH1D*)sparTrigger->Project(triPVZ);
        }
        else{
          sparTrigger->GetAxis(triMinv)->SetRangeUser(minfit,minRange-0.005);
          TH1D* low = (TH1D*)sparTrigger->Project(triPVZ);
          sparTrigger->GetAxis(triMinv)->SetRangeUser(maxRange+0.005,maxfit);
          TH1D* high = (TH1D*)sparTrigger->Project(triPVZ);
          hTriggersS = (TH1D*)low->Clone();
          hTriggersS->Add(high);
        }
        if(!hTriggersS) { printf("hTriggersS\n"); return; }
        hTriggersS->SetName(Form("TriggersS_%i_%i",ic,it));

        if(corr == "TPCFMDA"){
          sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(-5.5+0.001,-1.4-0.001);
          sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(-5.5+0.001,-1.4-0.001);
        }
        else{
          sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(1.4+0.001,4.0-0.001);
          sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(1.4+0.001,4.0-0.001);
        }
        TH2D* hPhiEtaMsum[3]={nullptr};
        TH2D* hPhiEtaSsum[3]={nullptr};
        Double_t nTriggersS =0.;
        Double_t nTriggersM =6666.;


        for (Int_t iz=1;iz<=nz;iz++){
          nTriggersS+=hTriggersS->Integral(iz,iz);
          sparSig->GetAxis(kVz)->SetRange(iz,iz);
          sparMix->GetAxis(kVz)->SetRange(iz,iz);

          Double_t zmin = sparSig->GetAxis(kVz)->GetBinLowEdge(iz);
          Double_t zmax = sparSig->GetAxis(kVz)->GetBinUpEdge(iz);

          TH2D *hPhiEtaS = nullptr;
          TH2D *hPhiEtaM = nullptr;
          if(iM == 0){
            sparSig->GetAxis(kMass)->SetRangeUser(minRange+0.005,maxRange-0.005);
            sparMix->GetAxis(kMass)->SetRangeUser(minRange+0.005,maxRange-0.005);
            hPhiEtaS  = (TH2D*)sparSig->Project(kdPhiTPCTPC,kdEtaTPCTPC);
            hPhiEtaM  = (TH2D*)sparMix->Project(kdPhiTPCTPC,kdEtaTPCTPC);
          }
          else{
            sparSig->GetAxis(kMass)->SetRangeUser(minfit,minRange-0.005);
            sparMix->GetAxis(kMass)->SetRangeUser(minfit,minRange-0.005);
            TH2D* lowSig = (TH2D*)sparSig->Project(kdPhiTPCTPC,kdEtaTPCTPC);
            TH2D* lowMix = (TH2D*)sparMix->Project(kdPhiTPCTPC,kdEtaTPCTPC);

            sparSig->GetAxis(kMass)->SetRangeUser(maxRange+0.005,maxfit);
            sparMix->GetAxis(kMass)->SetRangeUser(maxRange+0.005,maxfit);
            TH2D* highSig = (TH2D*)sparSig->Project(kdPhiTPCTPC,kdEtaTPCTPC);
            TH2D* highMix = (TH2D*)sparMix->Project(kdPhiTPCTPC,kdEtaTPCTPC);

            hPhiEtaS = (TH2D*)lowSig->Clone();
            hPhiEtaS->Add(highSig);
            hPhiEtaM = (TH2D*)lowMix->Clone();
            hPhiEtaM->Add(highMix);
          }
          if(!hPhiEtaS) { printf("hPhiEtaS\n"); return; }
          hPhiEtaS->SetName(Form("hPhiEtaS_%i_%i_%d_%i",ic,it,iM,iz));
          if(!hPhiEtaM) { printf("hPhiEtaM\n"); return; }
          hPhiEtaM->SetName(Form("hPhiEtaM_%i_%i_%d_%i",ic,it,iM,iz));

          // hPhiEtaS->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          // hPhiEtaM->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));


          Double_t norm = 1.;
          //====== For Norm
          Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-1*TMath::Pi()/2+0.0001);
          Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2-0.0001);

          // Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(0.-0.0001);
          // Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(0.+0.0001);

          Int_t binEta1, binEta2;

          if(corr == "TPCFMDA")
          {
            // -0.8 - -6. -> middle is -3.4
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(-2.6-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(-2.6+0.09);
          }
          else if(corr == "TPCFMDC"){
            // 1. - 4. -> middle is 2.5
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(2.5-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(2.5+0.09);
          }
          else if(corr == "FMDAFMDC"){
            // 3.4 - 8.2. -> middle is 5.8
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(5.8-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(5.8+0.09);
          }


          // Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-1.6-0.0001);
          // Int_t binEta1Ex = hPhiEtaM->GetYaxis()->FindBin(0.-0.01);
          // Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(1.6+0.0001);
          // Int_t binEta2Ex = hPhiEtaM->GetYaxis()->FindBin(0.0+0.0001);

          Double_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
          // Int_t nNormBins = (binEta2-binEta2Ex+1+binEta1Ex-binEta1+1)*(binPhi2-binPhi1+1);
          norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
          // norm = (hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta1Ex) + hPhiEtaM->Integral(binPhi1,binPhi2,binEta2Ex,binEta2))/nNormBins;
          if(norm < 1.0){
            printf("\n\n ******* \n\n");
            printf("bin eta : %d -- %d \n", binEta1, binEta2);
            printf("bin phi : %d -- %d \n", binPhi1, binPhi2);
            printf("bnorm : %f \n", norm);
          }

          //=======
          if(!hPhiEtaMsum[iM])
          {
           hPhiEtaMsum[iM] = (TH2D*) hPhiEtaM->Clone(Form("dphi_%i_%d_%i_M",it,iM,ic));
          }
          else
          {
           hPhiEtaMsum[iM]->Add(hPhiEtaM);
          }
          hPhiEtaM->Scale(1./norm);
          //hPhiEtaM->Write();
          // same/mixed
          TH2D* hPhiEtaSM = (TH2D*) hPhiEtaS->Clone(Form("hPhiEtaSM_%i_%i_%d_%i",ic,it,iM,iz));
          hPhiEtaSM->Divide(hPhiEtaM);
          if (!hPhiEtaSMsum[iM])
          {
           hPhiEtaSMsum[iM] = (TH2D*) hPhiEtaSM->Clone(Form("dphi_%i_%d_%i",it,iM,ic));
          }
          else
          {
           hPhiEtaSMsum[iM]->Add(hPhiEtaSM);
          }
          if (!hPhiEtaSsum[iM])
          {
           hPhiEtaSsum[iM] = (TH2D*)hPhiEtaS->Clone(Form("dphi_%i_%d_%i_S",it,iM,ic));
          }
          else
          {
           hPhiEtaSsum[iM]->Add(hPhiEtaS);
          }

        } //zvtx
        printf("nTriggersS=%f nTriggersM=%.0f\n",nTriggersS,nTriggersM);

        // if(useEff){
        //   Double_t efficiency = 0.;
        //   Int_t binEffStart = effHisto->GetXaxis()->FindBin(binst[it]+0.0001);
        //   Int_t binEffEnd = effHisto->GetXaxis()->FindBin(binst[it+1]-0.0001);
        //   for(Int_t i(binEffStart); i <= binEffEnd; i++){
        //     efficiency+= effHisto->GetBinContent(i);
        //   }
        //   efficiency = efficiency/(Double_t) (binEffEnd - binEffStart + 1);
        //   nTriggersS = nTriggersS/efficiency;
        //
        //   printf("scaling=%f nTriggersS=%.0f\n",efficiency,nTriggersS);
        // }


        nTrig[iM] = nTriggersS;

        hPhiEtaSMsum[iM]->Scale(1./nTriggersS);
        hPhiEtaSMsum[iM]->Scale(1./hPhiEtaSMsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum[iM]->Scale(1./hPhiEtaSMsum[iM]->GetYaxis()->GetBinWidth(1));

        hPhiEtaMsum[iM]->Scale(1./hPhiEtaMsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum[iM]->Scale(1./hPhiEtaMsum[iM]->GetYaxis()->GetBinWidth(1));

        // hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum[iM]->Scale(1./hPhiEtaSsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum[iM]->Scale(1./hPhiEtaSsum[iM]->GetYaxis()->GetBinWidth(1));
        hPhiEtaSMsum[iM]->SetTitle(Form("%.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,cen1,cen2));
        hPhiEtaSsum[iM]->SetTitle(Form("%.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,cen1,cen2));
        printf("tmax = %f\n",tmax);
        hPhiEtaSMsum[iM]->Write();
        hPhiEtaSsum[iM]->Write();
        hPhiEtaMsum[iM]->Write();



        // TCanvas* test = new TCanvas();
        // hPhiEtaSMsum->Draw("surf1");
        // test->SaveAs(Form("mix_%d.pdf",ia+1));
        //Draw Plots
          //       hPhiEtaSMsum->Draw("surf1 fb");return;
      } // iM end


      // Double_t signal = fitMinv->Integral(minRange+0.005,maxRange-0.005);
      // Double_t background = bkgMinv->Integral(minRange+0.005,maxRange-0.005);
      // nTrig[2] = (signal/(signal+background))*nTrig[0];
      hPhiEtaSMsum[2]=(TH2D*)hPhiEtaSMsum[0]->Clone();
      // hPhiEtaSMsum[2]->Add(hPhiEtaSMsum[1],-1);
      hPhiEtaSMsum[2]->SetName(Form("sub_%d",it));
      // hPhiEtaSMsum[2]->Scale(1./nTrig[2]);
      hPhiEtaSMsum[2]->Write();


    }   // pt -- trigger
    purity->Write();
 fout->Close();
 delete fout;

}
