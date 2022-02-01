#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

Double_t binst[]={1.8,4.8}; // FMDA
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

const Double_t binsa[]={-3.2,-1.8}; // FMDC
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

Double_t binscmin[] = {0.,60.};
Double_t binscmax[] = {20.,100.};
const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

Bool_t withAss = kFALSE;
Bool_t sameBin = kFALSE;
Bool_t useEff = kFALSE;
Bool_t isPP = kFALSE;

enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSamples};
enum {triPVZ, triSample};
// enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSample};

void readFMDFMD(Int_t sample = 0, TString corr = "FMDAFMDC", TString species = "Charged", Int_t ic = 0, TString cent = "peri", TString mcType = "gen")
// void readFMDFMD(TString corr = "FMDAFMDC", TString species = "Charged", Int_t ic = 0, TString cent = "peri", TString mcType = "gen")
// void readFMDFMD(TString corr = "1", TString species = "Charged", Int_t ic = 0, TString cent = "HM")
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

 const char* strPtTrg = "#eta_{FMD-A}";
 const char* strPtAss = "#eta_{FMD-C}";

 TFile *f = nullptr;
 // if(useEff) f = TFile::Open(Form("AR_eff_Nch_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 // else f = TFile::Open(Form("AR_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 // if(isPP) f = TFile::Open(Form("AR_pp_HM_%s.root",species.Data()),"READ");
 // f = TFile::Open(Form("tmp/AR_TPCFMDA_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_LHC16q/train_4229/AR_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data()),"READ");
 // f = TFile::Open(Form("/home/alidock/output/pp_LHC18/FMD_train_628_56/HM/AR_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data()),"READ");
 // f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_MC/train_1324/AR_%s_%s_%s.root",mcType.Data(),cent.Data(),species.Data()),"READ");
 if(!f) { printf("f\n"); return; }
 AliCFGridSparse *sparSig = (AliCFGridSparse*)f->Get("fhChargedSE_SelStep0");
 if(!sparSig) { printf("sparSig\n"); return; }

 // AliCFGridSparse *sparTriSig = (AliCFGridSparse*)f->Get("sparTriSig");
 AliCFGridSparse *sparMix = (AliCFGridSparse*)f->Get("fhChargedME_SelStep0");
 if(!sparMix) { printf("sparMix\n"); return; }
 // AliCFGridSparse *sparTriMix = (AliCFGridSparse*)f->Get("sparTriMix");
 // TH3D* trig = (TH3D*)f->Get("fhTrigTracks");
 // TH2D* trig = (TH2D*)f->Get("fhTrigTracks");
 // if(!trig) { printf("trig\n"); return; }
 // f->Close();

 AliCFGridSparse *sparTrigger = nullptr;
 if(species == "Charged") sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_SelStep0"));
 else sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));
 if(!sparTrigger) printf("Problem with sparTrigger! \n ");

 // AliCFGridSparse *sparTrigger = nullptr;
 // if(species == "Charged") sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_SelStep0"));
 // else sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));
 // if(!sparTrigger) printf("Problem with sparTrigger! \n ");

 // new TCanvas();
 // trig->Draw("colz");
 // return;

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = "";
 // if(useEff) outFileName = Form("Result_TPC_eff_%.0f_%.0f_%s_Nch_fullR.root",cen1,cen2,species.Data());
 // else outFileName = Form("Result_TPC_%.0f_%.0f_%s.root",cen1,cen2,species.Data());
 // if(isPP) outFileName = Form("Result_TPC_pp_HM_%s.root",species.Data());
 // outFileName = Form("Result_%s_%s_%s_testFMD.root",corr.Data(),cent.Data(),species.Data());
 // outFileName = Form("Result_%s_%s_%s_%s.root",mcType.Data(),corr.Data(),cent.Data(),species.Data());
 // outFileName = Form("Data/Samples_FMD/Result_%s_%s_%s_sample_%d.root",corr.Data(),cent.Data(),species.Data(),sample);
 outFileName = Form("Data/FMD_60_80/Result_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");


    for (Int_t it=0;it<nt;it++) {

      if(withAss && it == 1) continue;

      // sparSig->GetAxis(kSamples)->SetRangeUser(sample+0.001,sample+1.-0.001);
      // sparMix->GetAxis(kSamples)->SetRangeUser(sample+0.001,sample+1.-0.001);

      //pt or dphi trigger loop

      // sparTrigger->GetAxis(triSample)->SetRangeUser(sample+0.001,sample+1.-0.001);
      TH1D* hTriggersS = (TH1D*)sparTrigger->Project(triPVZ);
      // TH1D* hTriggersS = (TH1D*)trig->ProjectionY();
      // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d",it),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001));
      if(!hTriggersS) { printf("hTriggersS\n"); return; }
      hTriggersS->SetName(Form("TriggersS_%i_%i",ic,it));

      new TCanvas();
      hTriggersS->Draw();

      for (Int_t ia=0;ia<na;ia++) {
      //pt or eta assoc loop

        // if(corr == "1"){
        //   sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.6+0.001,5.1-0.001);
        //   sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.6+0.001,5.1-0.001);
        // }
        // if(corr == "2"){
        //   sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(4.4+0.001,6.-0.001);
        //   sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(4.4+0.001,6.-0.001);
        // }
        // if(corr == "3"){
        //   sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(5.7+0.001,6.5-0.001);
        //   sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(5.7+0.001,6.5-0.001);
        // }
        sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.6+0.001,6.8-0.001);
        sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.6+0.001,6.8-0.001);
        // // sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(4.6+0.001,8.-0.001);
        // sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(4.6+0.001,8.-0.001);
        // sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.7+0.001,8.-0.001);
        // sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(3.7+0.001,8.-0.001);

        // if(binst[it] < binsa[ia]) continue;

        if(sameBin && (binst[it] != binsa[ia])) continue;

        TH2D* hPhiEtaSMsum=0;
        TH2D* hPhiEtaSsum=0;
        TH2D* hPhiEtaMsum=0;
        Double_t nTriggersS =0.;
        Double_t nTriggersM =6666.;

        for (Int_t iz=1;iz<=nz;iz++){
          nTriggersS+=hTriggersS->Integral(iz,iz);
          sparSig->GetAxis(kVz)->SetRange(iz,iz);
          sparMix->GetAxis(kVz)->SetRange(iz,iz);

          Double_t zmin = sparSig->GetAxis(kVz)->GetBinLowEdge(iz);
          Double_t zmax = sparSig->GetAxis(kVz)->GetBinUpEdge(iz);

          TH2D *hPhiEtaS  = (TH2D*)sparSig->Project(kdPhiTPCTPC,kdEtaTPCTPC);
          if(!hPhiEtaS) { printf("hPhiEtaS\n"); return; }
          hPhiEtaS->SetName(Form("hPhiEtaS_%i_%i_%i_%i",ic,it,ia,iz));
          TH2D *hPhiEtaM  = (TH2D*)sparMix->Project(kdPhiTPCTPC,kdEtaTPCTPC);
          if(!hPhiEtaM) { printf("hPhiEtaM\n"); return; }
          hPhiEtaM->SetName(Form("hPhiEtaM_%i_%i_%i_%i",ic,it,ia,iz));

          // hPhiEtaS->Rebin2D(2.,2.);
          // hPhiEtaS->Scale(1./4.);
          // hPhiEtaM->Rebin2D(2.,2.);
          // hPhiEtaM->Scale(1./4.);

          hPhiEtaS->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          hPhiEtaM->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          //hPhiEtaS->Write();

          new TCanvas();
          hPhiEtaS->Draw("surf1");
          new TCanvas();
          hPhiEtaM->Draw("surf1");

          // if(iz < 3){
          //   TCanvas* can = new TCanvas();
          //   hPhiEtaS->Draw("surf1");
          //   can->SaveAs(Form("test_PVz_%d.pdf",iz));
          //
          //   TCanvas* can2 = new TCanvas();
          //   hPhiEtaM->Draw("surf1");
          //   can2->SaveAs(Form("test_PVz_ME_%d.pdf",iz));
          // }

          TCanvas* can = new TCanvas();
            hPhiEtaS->Draw("surf1");


          Double_t norm = 1.;
          //====== For Norm
          Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-1*TMath::Pi()/2+0.0001);
          Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2-0.0001);

          // Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(0.-0.0001);
          // Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(0.+0.0001);

          Int_t binEta1, binEta2;

          // if(corr == "TPCFMDA")
          // {
          //   // -0.8 - -6. -> middle is -3.4
          //   binEta1 = hPhiEtaM->GetYaxis()->FindBin(-3.4-0.09);
          //   binEta2 = hPhiEtaM->GetYaxis()->FindBin(-3.4+0.09);
          // }
          // else if(corr == "TPCFMDC"){
          //   // 1. - 4. -> middle is 2.5
          //   binEta1 = hPhiEtaM->GetYaxis()->FindBin(2.5-0.09);
          //   binEta2 = hPhiEtaM->GetYaxis()->FindBin(2.5+0.09);
          // }
          // else if(corr == "FMDAFMDC"){
          //   // 3.4 - 8.2. -> middle is 5.8
          //   binEta1 = hPhiEtaM->GetYaxis()->FindBin(6.0-0.09);
          //   binEta2 = hPhiEtaM->GetYaxis()->FindBin(6.0+0.09);
          // }
          if(corr == "1"){
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(4.1-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(4.1+0.09);
          }
          if(corr == "2"){
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(5.4-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(5.4+0.09);
          }
          if(corr == "3"){
            binEta1 = hPhiEtaM->GetYaxis()->FindBin(6.6-0.09);
            binEta2 = hPhiEtaM->GetYaxis()->FindBin(6.6+0.09);
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
            printf("it, iz, ia : %d, %d, %d \n", it, iz, ia);
          }

          //=======
          if(!hPhiEtaMsum)
          {
           hPhiEtaMsum = (TH2D*) hPhiEtaM->Clone(Form("dphi_%i_%i_%i_M",it,ia,ic));
          }
          else
          {
           hPhiEtaMsum->Add(hPhiEtaM);
          }

          hPhiEtaM->Scale(1./norm);
          //hPhiEtaM->Write();
          // same/mixed
          TH2D* hPhiEtaSM = (TH2D*) hPhiEtaS->Clone(Form("hPhiEtaSM_%i_%i_%i_%i",ic,it,ia,iz));
          hPhiEtaSM->Divide(hPhiEtaM);
          if (!hPhiEtaSMsum)
          {
           hPhiEtaSMsum = (TH2D*) hPhiEtaSM->Clone(Form("dphi_%i_%i_%i",it,ia,ic));
          }
          else
          {
           hPhiEtaSMsum->Add(hPhiEtaSM);
          }
          if (!hPhiEtaSsum)
          {
           hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphi_%i_%i_%i_S",it,ia,ic));
          }
          else
          {
           hPhiEtaSsum->Add(hPhiEtaS);
          }
          // new TCanvas();
          // hPhiEtaSM->Draw("surf1");

        } //zvtx
        printf("nTriggersS=%f nTriggersM=%.0f\n",nTriggersS,nTriggersM);

        hPhiEtaSMsum->Scale(1./nTriggersS);
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetYaxis()->GetBinWidth(1));

        // hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        // hPhiEtaSMsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        // hPhiEtaSsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        // printf("tmax = %f\n",tmax);
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();

        TCanvas* test = new TCanvas();
        hPhiEtaSMsum->Draw("surf1");
        // test->SaveAs(Form("mix_%d.pdf",ia+1));
        //Draw Plots
 //       hPhiEtaSMsum->Draw("surf1 fb");return;
      } // pt/eta -- associate
    }   // pt -- trigger
 // fout->Close();
 // delete fout;

}
