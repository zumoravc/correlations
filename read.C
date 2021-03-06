#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.5,0.6,0.7,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,10.0};
// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0};
// Double_t binst[]={0.2,0.5,1.0,1.5,2.0,3.0,4.0};
// Double_t binst[]={0.5,0.7,1.0,1.25,1.5,2.0,2.5,3.0,4.};
// Double_t binst[]={0.5, 1.0, 0.5, 0.7, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
// Double_t binst[]={0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
// Double_t binst[]={1.0, 1.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0};
// Double_t binst[]={1.5, 2.0};
Double_t binst[]={0.2, 3.0};
// Double_t binst[]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0};
//Double_t binst[]={0.5, 1.0, 1.5, 2.0, 5.0, 8.0};
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

// const Double_t binsa[]={0.5, 1., 1.5, 2.};
// const Double_t binsa[]={0.5, 1., 2., 200.};
// const Double_t binsa[]={1.5, 2.};
const Double_t binsa[]={0.2, 3.0};
// const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0};
// const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0};
// const Double_t binsa[]={0.2,0.5,1.0,1.5,2.0,3.0,4.0};
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

Double_t binscmin[] = {0.,60.};
Double_t binscmax[] = {20.,100.};
const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

Bool_t sameBin = kFALSE;
Bool_t useEff = kTRUE;
Bool_t isPP = kFALSE;

TString mcType = "reco";

TString fSystematicsFlagArr[5] = {"Ev0_Tr0", "Bayes", "Ev0_Tr1", "Ev2_Tr0", "Ev0_Tr9"};
TString sysTypeArr[5] = {"_default", "_Bayes", "_FB", "_PVz", "_TPC90"};

//correlation
enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSamples,kMass,kPt_TPC_trig};
// enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSamples,kPt_TPC_trig,kPt_ass};
//trigger
enum {triPVZ, triSample, triPt};

// void read(TString species = "Charged", Int_t flag = 3, Int_t ic = 1, Int_t sample = 0)
void read(TString species = "Charged", Int_t flag = 0, TString centName = "HM", Int_t ic = 0, Int_t sample = 0)
{
 TString fSystematicsFlag = fSystematicsFlagArr[flag];
 TString sysType = sysTypeArr[flag];

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

 // gStyle->SetOp tStat(0);
 const char* strPtTrg = "p_{T}^{Tri} (GeV/c)";
 const char* strPtAss = "p_{T}^{Asso} (GeV/c)";

 TFile *f = nullptr;
 // if(useEff) f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_LHC16q/train_4207/AR%s_%.0f_%.0f_%s.root",sysType.Data(),cen1,cen2,species.Data()),"READ");
 if(useEff) f = TFile::Open(Form("/Users/zuzana/Mirror/output/pp_merge/HM/dihcorr_1617/AR_%s_%s.root",centName.Data(),species.Data()),"READ");
 else f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_MC/train_1331/AR_%s_%.0f_%.0f_%s.root",mcType.Data(),cen1,cen2,species.Data()),"READ");
 if(!f) { printf("f\n"); return; }

 f->ls();

 AliCFGridSparse *sparSig = nullptr;
 if(species == "Charged") sparSig = (AliCFGridSparse*)f->Get("fhChargedSE_SelStep0");
 else sparSig = (AliCFGridSparse*)f->Get(Form("fhPidSE_%s_SelStep0",species.Data()));
 if(!sparSig) { printf("sparSig\n"); return; }

 AliCFGridSparse *sparMix = nullptr;
 if(species == "Charged") sparMix = (AliCFGridSparse*)f->Get("fhChargedME_SelStep0");
 else sparMix = (AliCFGridSparse*)f->Get(Form("fhPidME_%s_SelStep0",species.Data()));
 if(!sparMix) { printf("sparMix\n"); return; }

 AliCFGridSparse *sparTrigger = nullptr;
 if(species == "Charged") sparTrigger = (AliCFGridSparse*)f->Get("fhTrigTracks_SelStep0");
 else sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));

 // TH3D* trig = (TH3D*)f->Get(Form("fhTrigTracks_%s",species.Data()));

 //loading efficiencies
 TString part[4] = {"ch", "pi", "ka", "pr"};
 Int_t p = 0;
 if(species == "Pion") p = 1;
 else if(species == "Kaon") p = 2;
 else if(species == "Proton") p = 3;
 TH1D* effHisto = nullptr;
 TH2D* effHistoIn[21] = {nullptr};
 TH1D* effHistoInCent[21] = {nullptr};
 /*
 TFile *fileEff = TFile::Open("/Users/zuzana/alidock/correlations/Efficiencies/Efficiencies_pPb_syst.root", "READ");
 if(!fileEff) { printf("FileEff not open \n"); return; }
 TList* effList = (TList*)fileEff->Get("Efficiency2D_wFD");
 if(!effList) { printf("effList  not open \n"); return; }
 TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};
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
 */

 TFile *fileEff = TFile::Open("/Users/zuzana/alidock/correlations/Efficiencies/Efficiencies_pp_LHC16_HM_sys.root", "READ");
 TFile *fileEff2 = TFile::Open("/Users/zuzana/alidock/correlations/Efficiencies/Efficiencies_pp_LHC17_HM_sys.root", "READ");
 TList* effList = (TList*)fileEff->Get("Efficiency2D_wFD");
 // effList->ls();
 TList* effList2 = (TList*)fileEff2->Get("Efficiency2D_wFD");
 // effList2->ls();
 TString periods[10] = {"17e5", "17d3", "17f9", "17f5", "17d16", "17d17", "18f1", "17d18", "17f6", "18d8"};
 TString periods2[11] = {"18j4", "18a1", "17h1","17h11","17l5", "18c13", "18a9", "18a8", "18d3", "18c12", "17k4"};

 Int_t binCentMin = 0;
 Int_t binCentMax = 0;

 for(Int_t e(0); e < 10; e++){
   effHistoIn[e] = (TH2D*)effList->FindObject(Form("LHC%s_%s_%s_wFD",periods[e].Data(),"ch", fSystematicsFlag.Data()));
   if(!effHistoIn[e]) return;
   if(e == 0){
     binCentMin  = effHistoIn[e]->GetYaxis()->FindBin(0. + 0.001);
     binCentMax  = effHistoIn[e]->GetYaxis()->FindBin(100. - 0.001);
   }
   effHistoInCent[e] = (TH1D*) effHistoIn[e]->ProjectionX(Form("%s",part[p].Data()), binCentMin, binCentMax);
   effHistoInCent[e]->Scale(1./(binCentMax - binCentMin + 1));

   if(!effHisto) { effHisto = (TH1D*)effHistoInCent[e]->Clone("Efficiency"); effHisto->Sumw2(); }
   else effHisto->Add(effHistoInCent[e]);
  }
   for(Int_t e(0); e < 11; e++){
     effHistoIn[10+e] = (TH2D*)effList2->FindObject(Form("LHC%s_%s_%s_wFD",periods2[e].Data(),"ch", fSystematicsFlag.Data()));
     effHistoInCent[10+e] = (TH1D*) effHistoIn[e]->ProjectionX(Form("%s",part[p].Data()), binCentMin, binCentMax);
     effHistoInCent[10+e]->Scale(1./(binCentMax - binCentMin + 1));
     effHisto->Add(effHistoInCent[10+e]);
   }

  effHisto->Scale(1./21.);

  new TCanvas();
  effHisto->Draw();

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = "";
 // if(useEff) outFileName = Form("Data/Systematics/Result_TPC%s_eff_%.0f_%.0f_%s.root",sysType.Data(),cen1,cen2,species.Data());
 // if(useEff) outFileName = Form("Data/pp/Result_TPC_eff_%.0f_%.0f_%s_assRef_sample_%d.root",cen1,cen2,species.Data(),sample);
 if(useEff) outFileName = Form("Data/pp/Result_TPC_eff_%s_%s.root",centName.Data(),species.Data());
 else outFileName = Form("MC/Result_TPC_%s_%.0f_%.0f_%s_assRef.root",mcType.Data(),cen1,cen2,species.Data());
 // else outFileName = Form("Result_TPC_%.0f_%.0f_%s_noS_assRef.root",cen1,cen2,species.Data());
 if(isPP) outFileName = Form("pp/TPC_data/Result_TPC_pp_MB_%s_assRef.root",species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");


    for (Int_t it=0;it<nt;it++) {

      //pt or dphi trigger loop
      sparSig->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparMix->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);

      // sparSig->GetAxis(kSamples)->SetRangeUser(sample+0.001,sample+1.-0.001);
      // sparMix->GetAxis(kSamples)->SetRangeUser(sample+0.001,sample+1.-0.001);

      sparTrigger->GetAxis(triPt)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      // new TCanvas();
      // TH2D* triTest = (TH2D*)sparTrigger->Project(triPVZ,triSample);
      // triTest->Draw("colz");
      // return;
      // sparTrigger->GetAxis(triSample)->SetRangeUser(sample+0.001,sample+1.-0.001);

      Double_t tmin = sparSig->GetAxis(kPt_TPC_trig)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_trig)->GetFirst());
      Double_t tmax = sparSig->GetAxis(kPt_TPC_trig)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_trig)->GetLast());

      for (Int_t ia=0;ia<na;ia++) {
      //pt or eta assoc loop

        // if(binst[it] < binsa[ia]) continue;

        if(sameBin && (binst[it] != binsa[ia])) continue;

        TH2D* hPhiEtaSMsum=0;
        TH2D* hPhiEtaSsum=0;
        TH2D* hPhiEtaMsum=0;
        Double_t nTriggersS =0.;
        Double_t nTriggersM =6666.;
        for (Int_t iz=1;iz<=nz;iz++){
          // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d_sample_%d",it,iz),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001));
          TH1D* hTriggersS = (TH1D*)sparTrigger->Project(triPVZ);
          if(!hTriggersS) { printf("hTriggersS\n"); return; }
          hTriggersS->SetName(Form("TriggersS_%i_%i_%i",ic,it,iz));


          if(flag == 3){
            Double_t integ = hTriggersS->Integral(iz,iz);
            printf("Integral: %f \n", integ);
            if(integ < 1.) continue;
          }
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


          Double_t norm = 1.;
          //====== For Norm
          Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-1*TMath::Pi()/2+0.0001);
          Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2-0.0001);
          Int_t binPhi1Ex= hPhiEtaM->GetXaxis()->FindBin(-0.1);
          Int_t binPhi2Ex = hPhiEtaM->GetXaxis()->FindBin(0.1);

          // Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(0.-0.0001);
          // Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(0.+0.0001);

          Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(0.0-0.0001);
          Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(0.0+0.0001);
          //
          // Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-1.6-0.0001);
          // Int_t binEta1Ex = hPhiEtaM->GetYaxis()->FindBin(0.-0.01);
          // Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(1.6+0.0001);
          // Int_t binEta2Ex = hPhiEtaM->GetYaxis()->FindBin(0.0+0.0001);

          // Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
          // Int_t nNormBins = (binEta2-binEta2Ex+1+binEta1Ex-binEta1+1)*(binPhi2-binPhi1+1);
          Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi2Ex+1+binPhi1Ex-binPhi1+1);
          // norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
          // norm = (hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta1Ex) + hPhiEtaM->Integral(binPhi1,binPhi2,binEta2Ex,binEta2))/nNormBins;
          norm = (hPhiEtaM->Integral(binPhi1,binPhi1Ex,binEta1,binEta2) + hPhiEtaM->Integral(binPhi2Ex,binPhi2,binEta1,binEta2))/nNormBins;
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

        // hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSMsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        hPhiEtaSsum->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        printf("tmax = %f\n",tmax);
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
