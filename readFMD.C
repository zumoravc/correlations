#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
Double_t binst[]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
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
Bool_t useEff = kTRUE;
Bool_t isPP = kFALSE;

TString fSystematicsFlagArr[5] = {"Ev0_Tr0", "Bayes", "Ev0_Tr1", "Ev2_Tr0", "Ev0_Tr9"};
TString sysTypeArr[5] = {"_default", "_Bayes", "_FB", "_PVz", "_TPC90"};


enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSample,kMass,kPt_TPC_trig};
enum {triPVZ, triSample, triPt};

void readFMD(Int_t sample = 0, TString corr = "TPCFMDA", Int_t flag = 0, TString species = "Kaon", Int_t ic = 0, TString cent = "cent", TString mcType = "gen")
// void readFMD(TString corr = "TPCFMDA", Int_t flag = 0, TString species = "Pion", Int_t ic = 0, TString cent = "cent", TString mcType = "gen")
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

 TString fSystematicsFlag = fSystematicsFlagArr[flag];
 TString sysType = sysTypeArr[flag];


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

 TFile *f = nullptr;
 // if(useEff) f = TFile::Open(Form("AR_eff_Nch_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 // else f = TFile::Open(Form("AR_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 // if(isPP) f = TFile::Open(Form("AR_pp_HM_%s.root",species.Data()),"READ");
 // f = TFile::Open(Form("tmp/AR_TPCFMDA_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_LHC16q/train_4228/AR_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data()),"READ");
 // f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_MC/train_1328/AR_%s_%s_%s.root",mcType.Data(),cent.Data(),species.Data()),"READ");
 // f = TFile::Open(Form("/home/alidock/output/pp_LHC18/train_6247/test/AR_%s_%s_%s.root",corr.Data(),cent.Data(),species.Data()),"READ");
 if(!f) { printf("f\n"); return; }
 AliCFGridSparse *sparSig = nullptr;
 if(species == "Charged") sparSig = (AliCFGridSparse*)f->Get("fhChargedSE_SelStep0");
 else if(species == "Pion" || species == "Kaon" || species == "Proton") sparSig = (AliCFGridSparse*)f->Get(Form("fhPidSE_%s_SelStep0",species.Data()));
 else printf("***** This should not happen. ***** \n");
 if(!sparSig) { printf("sparSig\n"); return; }

 // AliCFGridSparse *sparTriSig = (AliCFGridSparse*)f->Get("sparTriSig");
 AliCFGridSparse *sparMix = nullptr;
 if(species == "Charged") sparMix = (AliCFGridSparse*)f->Get("fhChargedME_SelStep0");
 else if(species == "Pion" || species == "Kaon" || species == "Proton") sparMix = (AliCFGridSparse*)f->Get(Form("fhPidME_%s_SelStep0",species.Data()));
 else printf("***** This should not happen. ***** \n");
 if(!sparMix) { printf("sparMix\n"); return; }
 // AliCFGridSparse *sparTriMix = (AliCFGridSparse*)f->Get("sparTriMix");
 // TH3D* trig = nullptr;
 // if(species == "Charged") trig = (TH3D*)f->Get("fhTrigTracks");
 // else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH3D*)f->Get(Form("fhTrigTracks_%s",species.Data()));
 // if(!trig) { printf("trig\n"); return; }
 // TH2D* trig = nullptr;
 // if(species == "Charged") trig = (TH2D*)f->Get("fhTrigTracks");
 // else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH2D*)f->Get(Form("fhTrigTracks_%s",species.Data()));
 // if(!trig) { printf("trig\n"); return; }
 // f->Close();

 AliCFGridSparse *sparTrigger = nullptr;
 if(species == "Charged") sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_SelStep0"));
 else sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));
 if(!sparTrigger) printf("Problem with sparTrigger! \n ");

 // TFile *fileEff = TFile::Open("task/Efficiency_LHC20f11c_NoFD_PID.root", "READ");
 // if(!fileEff) { printf("FileEff  not open \n"); return; }
 // TList* effList = (TList*)fileEff->Get("EffAndFD");
 // if(!effList) { printf("effList  not open \n"); return; }
 // TH1D* effHisto = nullptr;
 // if(species == "Charged") effHisto = (TH1D*)effList->FindObject("EffRescaled_ch_eta0");
 // else if(species == "Pion") effHisto = (TH1D*)effList->FindObject("EffRescaled_pi_eta0");
 // else if(species == "Kaon") effHisto = (TH1D*)effList->FindObject("EffRescaled_ka_eta0");
 // else effHisto = (TH1D*)effList->FindObject("EffRescaled_pr_eta0");
 // if(!effHisto) { printf("effHisto not open \n"); return; }

 TFile *fileEff = TFile::Open("/Users/zuzana/alidock/correlations/Efficiencies/Efficiencies_pPb_syst.root", "READ");
 if(!fileEff) { printf("FileEff not open \n"); return; }
 TList* effList = (TList*)fileEff->Get("Efficiency2D_wFD");
 if(!effList) { printf("effList  not open \n"); return; }
 TH2D* effHistoIn[8] = {nullptr};
 TH1D* effHistoInCent[8] = {nullptr};
 TH1D* effHisto = nullptr;
 TString part[4] = {"ch", "pi", "ka", "pr"};
 TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
 Int_t p = 0;
 if(species == "Pion") p = 1;
 else if(species == "Kaon") p = 2;
 else if(species == "Proton") p = 3;
 for(Int_t e(0); e < 8; e++){
   effHistoIn[e] = (TH2D*)effList->FindObject(Form("LHC17f2b_%s_Eta_%s_%s_wFD",part[p].Data(), etaReg[e].Data(),fSystematicsFlag.Data()));
   if(!effHistoIn[e]) { printf("Efficiency (%s, eta region %s, flag %s) not loaded",part[p].Data(),etaReg[e].Data(),fSystematicsFlag.Data()); return; }

   Int_t binCentMin  = effHistoIn[e]->GetYaxis()->FindBin(cen1 + 0.001);
   Int_t binCentMax  = effHistoIn[e]->GetYaxis()->FindBin(cen2 - 0.001);

   effHistoInCent[e] = (TH1D*) effHistoIn[e]->ProjectionX(Form("%s_Eta_%s",part[p].Data(), etaReg[e].Data()), binCentMin, binCentMax);
   effHistoInCent[e]->Scale(1./(binCentMax - binCentMin + 1));

   if(!effHisto) { effHisto = (TH1D*)effHistoInCent[e]->Clone("Efficiency"); effHisto->Sumw2(); }
   else effHisto->Add(effHistoInCent[e]);
 }
 effHisto->Scale(1./8.);

 // new TCanvas();
 // trig->Draw("colz");
 // return;

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = "";
 if(useEff) outFileName = Form("Data/Samples_FMD/Result_%s_eff_%s_%s_sample_%d.root",corr.Data(),cent.Data(),species.Data(),sample);
 // if(useEff) outFileName = Form("Data/Result_%s_eff_%s_%s.root",corr.Data(),cent.Data(),species.Data());
 else outFileName = Form("Result_%s_%s_%s_%s.root",mcType.Data(),corr.Data(),cent.Data(),species.Data());
 // if(isPP) outFileName = Form("Result_TPC_pp_HM_%s.root",species.Data());
 // outFileName = Form("Result_%s_%s_%s_new.root",corr.Data(),cent.Data(),species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");


    for (Int_t it=0;it<nt;it++) {

      if(withAss && it == 1) continue;

      //pt or dphi trigger loop
      sparSig->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparMix->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);

      sparSig->GetAxis(kSample)->SetRangeUser(sample+0.001,sample+1.-0.001);
      sparMix->GetAxis(kSample)->SetRangeUser(sample+0.001,sample+1.-0.001);

      // sparTriSig->GetAxis(kPt_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      //sparTriMix->GetAxis(kPt_trig_Mix)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      Double_t tmin = sparSig->GetAxis(kPt_TPC_trig)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_trig)->GetFirst());
      Double_t tmax = sparSig->GetAxis(kPt_TPC_trig)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_trig)->GetLast());

      sparTrigger->GetAxis(triPt)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparTrigger->GetAxis(triSample)->SetRangeUser(sample+0.001,sample+1.-0.001);
      TH1D* hTriggersS = (TH1D*)sparTrigger->Project(triPVZ);
      // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d",it),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001));
      if(!hTriggersS) { printf("hTriggersS\n"); return; }
      hTriggersS->SetName(Form("TriggersS_%i_%i",ic,it));

      for (Int_t ia=0;ia<na;ia++) {
      //pt or eta assoc loop

        // if(binst[it] < binsa[ia]) continue;

        if(sameBin && (binst[it] != binsa[ia])) continue;

        if(corr == "TPCFMDA"){
          sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(-5.5+0.001,-1.4-0.001);
          sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(-5.5+0.001,-1.4-0.001);
        }
        else{
          sparSig->GetAxis(kdEtaTPCTPC)->SetRangeUser(1.4+0.001,4.0-0.001);
          sparMix->GetAxis(kdEtaTPCTPC)->SetRangeUser(1.4+0.001,4.0-0.001);
        }
        // Double_t amin = sparSig->GetAxis(kPt_TPC_asso)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_asso)->GetFirst());
        // Double_t amax = sparSig->GetAxis(kPt_TPC_asso)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_asso)->GetLast());

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

          hPhiEtaS->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          hPhiEtaM->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          //hPhiEtaS->Write();

          // hPhiEtaS->Rebin2D(2.,2.);
          // hPhiEtaS->Scale(1./4.);
          // hPhiEtaM->Rebin2D(2.,2.);
          // hPhiEtaM->Scale(1./4.);

          // if(iz < 3){
          //   TCanvas* can = new TCanvas();
          //   hPhiEtaS->Draw("surf1");
          //   can->SaveAs(Form("test_PVz_%d.pdf",iz));
          //
          //   TCanvas* can2 = new TCanvas();
          //   hPhiEtaM->Draw("surf1");
          //   can2->SaveAs(Form("test_PVz_ME_%d.pdf",iz));
          // }


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

        } //zvtx
        printf("nTriggersS=%f nTriggersM=%.0f\n",nTriggersS,nTriggersM);

        if(useEff){
          Double_t efficiency = 0.;
          Int_t binEffStart = effHisto->GetXaxis()->FindBin(binst[it]+0.0001);
          Int_t binEffEnd = effHisto->GetXaxis()->FindBin(binst[it+1]-0.0001);
          for(Int_t i(binEffStart); i <= binEffEnd; i++){
            efficiency+= effHisto->GetBinContent(i);
          }
          efficiency = efficiency/(Double_t) (binEffEnd - binEffStart + 1);
          nTriggersS = nTriggersS/efficiency;

          printf("scaling=%f nTriggersS=%.0f\n",efficiency,nTriggersS);
        }


        hPhiEtaSMsum->Scale(1./nTriggersS);
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSMsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        hPhiEtaSsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        printf("tmax = %f\n",tmax);
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();

        // TCanvas* test = new TCanvas();
        // hPhiEtaSMsum->Draw("surf1");
        // test->SaveAs(Form("mix_%d.pdf",ia+1));
        //Draw Plots
 //       hPhiEtaSMsum->Draw("surf1 fb");return;
      } // pt/eta -- associate
    }   // pt -- trigger
 fout->Close();
 delete fout;

}
