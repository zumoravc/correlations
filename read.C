#include "AliTHn.h"
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";

Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binst[]={0.5, 1.0, 0.5, 0.7, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
// Double_t binst[]={0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
// Double_t binst[]={1.0, 1.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0};
// Double_t binst[]={1.0, 2.0, 1.0, 2.0, 4.0};
//Double_t binst[]={0.5, 1.0, 1.5, 2.0, 5.0, 8.0};
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

// const Double_t binsa[]={0.5, 1., 1.5, 2.};
// const Double_t binsa[]={0.5, 1., 2., 200.};
// const Double_t binsa[]={1., 2.};
// const Double_t binsa[]={0.5, 1.0};
const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

Double_t binscmin[] = {0.,60.};
Double_t binscmax[] = {20.,100.1};
const Int_t nc = sizeof(binscmin) / sizeof(Double_t);

Bool_t withAss = kFALSE;
Bool_t sameBin = kTRUE;

//signal correlation
// enum {kdEtaTPCTPC,kPt_TPC_asso,kPt_TPC_trig,kdPhiTPCTPC,kVz};
enum {kdEtaTPCTPC,kdPhiTPCTPC,kPt_TPC_trig,kPt_TPC_asso,kVz};
//signal trigger
// enum {kPt_trig,kVz_Tri,kRan};
//mix trigger
//enum {kPt_trig_Mix,kCen_Tri_Mix};

void read(TString species = "Pion", Int_t ic = 0)
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

 TFile *f = TFile::Open(Form("AR_%.0f_%.0f_%s.root",cen1,cen2,species.Data()),"READ");
 // TFile *f = TFile::Open("corr_p.root","READ");
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
 TH2D* trig = nullptr;
 if(species == "Charged") trig = (TH2D*)f->Get("fhTrigTracks");
 else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH2D*)f->Get(Form("fhTrigTracks_%s",species.Data()));
 if(!trig) { printf("trig\n"); return; }
 // f->Close();

 // new TCanvas();
 // trig->Draw("colz");
 // return;

 Int_t nz = sparSig->GetAxis(kVz)->GetNbins();
 printf("zbins = %i\n",nz);

 TString outFileName = Form("Result_TPC_%.0f_%.0f_%s.root",cen1,cen2,species.Data());
 TFile* fout = TFile::Open(outFileName,"RECREATE");


    for (Int_t it=0;it<nt;it++) {

      if(withAss && it == 1) continue;

      //pt or dphi trigger loop
      sparSig->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      sparMix->GetAxis(kPt_TPC_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      // sparTriSig->GetAxis(kPt_trig)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      //sparTriMix->GetAxis(kPt_trig_Mix)->SetRangeUser(binst[it]+0.001,binst[it+1]-0.001);
      Double_t tmin = sparSig->GetAxis(kPt_TPC_trig)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_trig)->GetFirst());
      Double_t tmax = sparSig->GetAxis(kPt_TPC_trig)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_trig)->GetLast());

      TH1D* hTriggersS = (TH1D*)trig->ProjectionY(Form("py_%d",it),trig->GetXaxis()->FindBin(binst[it]+0.001),trig->GetXaxis()->FindBin(binst[it+1]-0.001));
      if(!hTriggersS) { printf("hTriggersS\n"); return; }
      hTriggersS->SetName(Form("TriggersS_%i_%i",ic,it));

      for (Int_t ia=0;ia<na;ia++) {
      //pt or eta assoc loop

        if(binst[it] < binsa[ia]) continue;

        if(sameBin && (binst[it] != binsa[ia])) continue;

        sparSig->GetAxis(kPt_TPC_asso)->SetRangeUser(binsa[ia]+0.001,binsa[ia+1]-0.001);
        sparMix->GetAxis(kPt_TPC_asso)->SetRangeUser(binsa[ia]+0.001,binsa[ia+1]-0.001);
        Double_t amin = sparSig->GetAxis(kPt_TPC_asso)->GetBinLowEdge(sparSig->GetAxis(kPt_TPC_asso)->GetFirst());
        Double_t amax = sparSig->GetAxis(kPt_TPC_asso)->GetBinUpEdge (sparSig->GetAxis(kPt_TPC_asso)->GetLast());

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

          // TCanvas* test = new TCanvas();
          // hPhiEtaM->Draw("surf1");
          // test->SaveAs(Form("test_PVz_%d.pdf",iz));

          Double_t norm = 1.;
          //====== For Norm
          // Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-1*TMath::Pi()/2+0.0001);
          // Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2-0.0001);

          Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(0.-0.0001);
          Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(0.+0.0001);

          Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(0.-0.0001);
          Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(0.+0.0001);

          Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
          norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
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
        printf("nTriggersS=%.0f nTriggersM=%.0f\n",nTriggersS,nTriggersM);
        hPhiEtaSMsum->Scale(1./nTriggersS);
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1./hPhiEtaMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaSsum->Scale(1./nTriggersS);
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1./hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSMsum->SetTitle(Form("%.2f < %s < %.2f , %.1f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
        hPhiEtaSsum->SetTitle(Form("%.2f < %s < %.2f , %.1f < %s < %.2f , %.0f~%.0f%%",tmin,strPtTrg,tmax,binsa[ia],strPtAss,binsa[ia+1],cen1,cen2));
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
