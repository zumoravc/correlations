#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
Double_t binst[]={0.5,0.7,0.9,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.6,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={1.25,1.5,1.75};
// Double_t binst[]={0.5,0.6};
// Double_t binst[]={0.2, 3.0};
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

const Double_t binsa[]={0.2, 3.0};
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

Double_t binscmin[] = {0.,60.};
Double_t binscmax[] = {20.,100.};
const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

Bool_t sameBin = kFALSE;
Bool_t useEff = kFALSE;
Bool_t isPP = kFALSE;

//signal correlation
// enum {kdEtaTPCTPC,kPt_TPC_asso,kPt_TPC_trig,kdPhiTPCTPC,kVz};
enum {kdEtaTPCTPC,kdPhiTPCTPC,kVz,kSamples,kMass,kPt_TPC_trig,kPt_TPC_asso};
// enum {kdEtaTPCTPC,kdPhiTPCTPC,kPt_TPC_trig,kPt_TPC_asso,kVz};
//signal trigger
// enum {kPt_trig,kVz_Tri,kRan};
//mix trigger
//enum {kPt_trig_Mix,kCen_Tri_Mix};
enum {triMinv, triPt, triPVZ, triSample};

void readV0(TString species = "K0s", Int_t ic = 0, Int_t sample = 0)
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
 //Double_t binscmin[] = {0., 0,  0., 60., 60., 0.};
 Double_t binscmax[] = {cen2};
 //Double_t binscmax[] = {5., 10, 20., 70., 100.1, 100.1};
 const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

 // gStyle->SetOp tStat(0);
 const char* strPtTrg = "p_{T}^{Tri} (GeV/c)";
 const char* strPtAss = "p_{T}^{Asso} (GeV/c)";

 TFile *f = nullptr;
 if(useEff) f = TFile::Open(Form("AR_eff_Nch_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 else f = TFile::Open(Form("/Users/zuzana/Mirror/output/pPb_LHC16q/train_4189/AR__%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 if(isPP) f = TFile::Open(Form("AR_TPC_pp_MB_%s.root",species.Data()),"READ");
 // TFile *f = TFile::Open("corr_p.root","READ");
 if(!f) { printf("f\n"); return; }
 AliCFGridSparse *sparSig = nullptr;
 if(species == "Charged") sparSig = (AliCFGridSparse*)f->Get("fhChargedSE_SelStep0");
 else if(species == "Pion" || species == "Kaon" || species == "Proton" || species == "K0s" || species == "Lambda") sparSig = (AliCFGridSparse*)f->Get(Form("fhPidSE_%s_SelStep0",species.Data()));
 else printf("***** This should not happen. ***** \n");
 if(!sparSig) { printf("sparSig\n"); return; }

 // AliCFGridSparse *sparTriSig = (AliCFGridSparse*)f->Get("sparTriSig");
 AliCFGridSparse *sparMix = nullptr;
 if(species == "Charged") sparMix = (AliCFGridSparse*)f->Get("fhChargedME_SelStep0");
 else if(species == "Pion" || species == "Kaon" || species == "Proton" || species == "K0s" || species == "Lambda") sparMix = (AliCFGridSparse*)f->Get(Form("fhPidME_%s_SelStep0",species.Data()));
 else printf("***** This should not happen. ***** \n");
 if(!sparMix) { printf("sparMix\n"); return; }
 // AliCFGridSparse *sparTriMix = (AliCFGridSparse*)f->Get("sparTriMix");
 // TH2D* trig = nullptr;
 // if(species == "Charged") trig = (TH2D*)f->Get("fhTrigTracks");
 // else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH2D*)f->Get(Form("fhTrigTracks_%s",species.Data()));
 TH3D* trig = nullptr;
 if(species == "Charged") trig = (TH3D*)f->Get("fhTrigTracks");
 else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH3D*)f->Get(Form("fhTrigTracks_%s",species.Data()));
 // if(!trig) { printf("trig\n"); return; }
 // f->Close();


 AliCFGridSparse *sparTrigger = nullptr;
 if(species == "K0s" || species == "Lambda") sparTrigger = (AliCFGridSparse*)f->Get(Form("fhTrigTracks_%s_SelStep0",species.Data()));
 if(!sparTrigger) printf("Problem with sparTrigger! \n ");

 TFile *fileEff = TFile::Open("task/Efficiency_LHC20f11c_NoFD_PID.root", "READ");
 if(!fileEff) { printf("FileEff  not open \n"); return; }
 TList* effList = (TList*)fileEff->Get("EffAndFD");
 if(!effList) { printf("effList  not open \n"); return; }
 TH1D* effHisto = nullptr;
 if(species == "Charged") effHisto = (TH1D*)effList->FindObject("EffRescaled_ch_eta0");
 else if(species == "Pion") effHisto = (TH1D*)effList->FindObject("EffRescaled_pi_eta0");
 else if(species == "Kaon") effHisto = (TH1D*)effList->FindObject("EffRescaled_ka_eta0");
 else effHisto = (TH1D*)effList->FindObject("EffRescaled_pr_eta0");
 if(!effHisto) { printf("effHisto not open \n"); return; }

 Int_t nTrig[3] = {0,0,0};

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = "";
 if(useEff) outFileName = Form("Result_TPC_eff_%.0f_%.0f_%s_Nch_fullR.root",cen1,cen2,species.Data());
 else outFileName = Form("Result_TPC_%.0f_%.0f_%s_assRef_new.root",cen1,cen2,species.Data());
 // else outFileName = Form("Result_TPC_%.0f_%.0f_%s_noS_assRef.root",cen1,cen2,species.Data());
 if(isPP) outFileName = Form("pp/TPC_data/Result_TPC_pp_MB_%s_assRef.root",species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");

 if(species != "K0s" && species != "Lambda") {printf("NOT V0!!! \n"); return; }

    for (Int_t it=0;it<nt;it++) {

      Double_t minfit = 0.0;
      Double_t maxfit = 0.0;
      if(species == "K0s"){
        minfit = 0.44;
        maxfit = 0.56;
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



      // printf("mean %f \t %f \n", fitMinv->GetParameter(6), fitMinv->GetParameter(8));
      // printf("sigma %f \t %f \n", fitMinv->GetParameter(7), fitMinv->GetParameter(9));

      Double_t meanPDG = 0.0;
      if(species == "K0s") meanPDG = 0.497614;
      else meanPDG = 1.11568;

      Double_t sigmaFit = fitMinv->GetParameter(7);
      if(fitMinv->GetParameter(9) > sigmaFit) sigmaFit = fitMinv->GetParameter(9);

      Double_t minRange = meanPDG - 3.0*sigmaFit;
      Double_t maxRange = meanPDG + 3.0*sigmaFit;

      printf("\n range: %f -- %f \n\n", minRange, maxRange);


      //pt or dphi trigger loop
      sparSig->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparMix->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);

      Double_t tmin = sparSig->GetAxis(kPt_TPC_trig)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_trig)->GetFirst());
      Double_t tmax = sparSig->GetAxis(kPt_TPC_trig)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_trig)->GetLast());

      for (Int_t ia=0;ia<na;ia++) {

        // if(binst[it] < binsa[ia]) continue;

        if(sameBin && (binst[it] != binsa[ia])) continue;

        sparSig->GetAxis(kPt_TPC_asso)->SetRangeUser(binsa[ia]+0.001,binsa[ia+1]-0.001);
        sparMix->GetAxis(kPt_TPC_asso)->SetRangeUser(binsa[ia]+0.001,binsa[ia+1]-0.001);

        TH2D* hPhiEtaSMsum[3]={nullptr};
        TH2D* hPhiEtaSsum[3]={nullptr};
        TH2D* hPhiEtaMsum[3]={nullptr};
        Double_t nTriggersS =0.;
        Double_t nTriggersM =6666.;

        for(Int_t iM(0); iM < 2; iM++){
          if(iM == 1) continue;

        for (Int_t iz=1;iz<=nz;iz++){
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

          // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d_sample_%d",it,iz),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001));
          // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d_sample_%d",it,iz),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001),trig->GetZaxis()->FindBin(sample+0.001),trig->GetZaxis()->FindBin(sample+1-0.001));
          // TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d_sample_%d",it,iz),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001),trig->GetZaxis()->FindBin(sample+0.001),trig->GetZaxis()->FindBin(9+1-0.001));
          if(!hTriggersS) { printf("hTriggersS\n"); return; }
          hTriggersS->SetName(Form("TriggersS_%i_%i_%i",ic,it,iz));

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

          hPhiEtaS->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          hPhiEtaM->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",binsa[ia],strPtAss,binsa[na],zmin,zmax));
          //hPhiEtaS->Write();

          // TCanvas* test = new TCanvas();
          // hPhiEtaM->Draw("surf1");
          // test->SaveAs(Form("test_PVz_%d.pdf",iz));

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

        nTrig[iM] = nTriggersS;

        hPhiEtaSMsum[iM]->Scale(1./nTriggersS);
        hPhiEtaSMsum[iM]->Scale(1./hPhiEtaSMsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum[iM]->Scale(1./hPhiEtaSMsum[iM]->GetYaxis()->GetBinWidth(1));

        hPhiEtaMsum[iM]->Scale(1./hPhiEtaMsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum[iM]->Scale(1./hPhiEtaMsum[iM]->GetYaxis()->GetBinWidth(1));

        // hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum[iM]->Scale(1./hPhiEtaSsum[iM]->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum[iM]->Scale(1./hPhiEtaSsum[iM]->GetYaxis()->GetBinWidth(1));
        hPhiEtaSMsum[iM]->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        hPhiEtaSsum[iM]->SetTitle(Form("%.2f < %s < %.2f , %.2f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        printf("tmax = %f\n",tmax);
        hPhiEtaSMsum[iM]->Write();
        hPhiEtaSsum[iM]->Write();
        hPhiEtaMsum[iM]->Write();

        TCanvas* test = new TCanvas();
        hPhiEtaSMsum[iM]->Draw("surf1");
      } // end iM

      // Double_t signal = fitMinv->Integral(minRange+0.005,maxRange-0.005);
      // Double_t background = bkgMinv->Integral(minRange+0.005,maxRange-0.005);
      // nTrig[2] = (signal/(signal+background))*nTrig[0];
      hPhiEtaSMsum[2]=(TH2D*)hPhiEtaSMsum[0]->Clone();
      // hPhiEtaSMsum[2]->Add(hPhiEtaSMsum[1],-1);
      hPhiEtaSMsum[2]->SetName(Form("sub_%d",it));
      // hPhiEtaSMsum[2]->Scale(1./nTrig[2]);
      hPhiEtaSMsum[2]->Write();

        // test->SaveAs(Form("mix_%d.pdf",ia+1));
        //Draw Plots
 //       hPhiEtaSMsum->Draw("surf1 fb");return;
      } // pt/eta -- associate
    }   // pt -- trigger
 // fout->Close();
 // delete fout;

}
