#include <TPad.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include "./env.h"


const char* kTitlePhi          = "#Delta#varphi (rad)";
const char* kTitleEta          = "#Delta#eta";
const char* kCorrFuncTitle =       "Y_{sub} (rad^{-1})";
const char* kCorrFuncTitle_Before =       "Y (rad^{-1})";
const char* kProjYieldTitlePhi = "Y_{sub} per #Delta#eta (rad^{-1})";
const Double_t pi = TMath::Pi();
void SymmetrizeTo0toPi(TH1D* h);
void PadFor2DCorr();
TH1D* LinearFitProjectionX(TH2D* h2, TString name, const Int_t bin1, const Int_t bin2, Bool_t UseStandardProj= kFALSE,Int_t debug=0);
Double_t fun_template(Double_t *x, Double_t *par);
void minuitfcn23_1(int &npar, Double_t *gin, Double_t &ff, Double_t *par,int iflag);
void tempminuit23_1(Double_t*fParamVal,Double_t*fParamErr);

TH1D* periH;

//for template fit
TH1D* hminuit_peri[3];
TH1D* hminuit[3];

Int_t counterCorr = 0;

// Double_t binst[]={0.2,0.3,0.5,0.7,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
// Double_t binst[]={0.2,3.0};
// Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
Double_t binst[]={0.5,0.7,0.9,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
// Double_t binst[]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

const Double_t binsa[]={0.2,3.0};
// const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

// Double_t binsV2[11] = {0.2,0.3,0.5,0.7,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
// Double_t binsV2[22] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binsV2[2] = {0.2,3.};
// Double_t binsV2[16] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
Double_t binsV2[14] = {0.5,0.7,0.9,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};

Bool_t withAss = kFALSE;
Bool_t sameBin = kFALSE;
Bool_t useCentReg = kFALSE;
Bool_t useZYAM = kFALSE;
Bool_t useAssRef = kFALSE;
Bool_t combine = kTRUE;

void templatefitFMD(        TString species = "Lambda",
                  Int_t sample = 0,
                  // TString folder = "plotsWithEfficiency/pPb/FMD",
                  TString folder = "pPb/FMD_V0_test",
                  TString sStyle = "CC",
                  TString sEst = "V0A",
                  TString mcType = "reco",
                  Int_t c1 = 0,
                  Int_t c2 = 1,
                  Double_t* res = NULL)
{

  TString corr[9] = {"TPCFMDA", "TPCFMDC", "FMDAFMDC", "TPCFMDA", "TPCFMDC", "FMDAFMDC", "TPCFMDA", "TPCFMDC", "FMDAFMDC"};
  TString centNames[9] = {"cent", "cent", "cent", "peri", "peri", "peri", "sub", "sub", "sub"};
  TString centLegend[9] = {"0-20%", "0-20%", "0-20%", "60-100%", "60-100%", "60-100%", "(0-20%)-(60-100%)", "(0-20%)-(60-100%)", "(0-20%)-(60-100%)"};
  // TString centNames[9] = {"HM", "HM", "HM", "MB", "MB", "MB", "sub", "sub", "sub"};
  // TString centLegend[9] = {"HM", "HM", "HM", "MB", "MB", "MB", "HM-MB", "HM-MB", "HM-MB"};


  TFile* file[6]; // first three from HM, then from LM
  for(Int_t i(0); i < 6; i++){
    // if(i != 2 && i != 5)  file[i] = TFile::Open(Form("Data/Result_%s_eff_%s_%s.root",corr[i].Data(),centNames[i].Data(),species.Data()),"READ");
    if(i != 2 && i != 5)  file[i] = TFile::Open(Form("FMD_Data/Result_%s_%s_%s_test.root",corr[i].Data(),centNames[i].Data(),species.Data()),"READ");
    else file[i] = TFile::Open(Form("Data/Result_%s_%s_%s.root",corr[i].Data(),centNames[i].Data(),"Charged"),"READ");
    // if(i != 2 && i != 5)  file[i] = TFile::Open(Form("Data/Samples_FMD/Result_%s_eff_%s_%s_sample_%d.root",corr[i].Data(),centNames[i].Data(),species.Data(),sample),"READ");
    // else file[i] = TFile::Open(Form("Data/Samples_FMD/Result_%s_%s_%s_sample_%d.root",corr[i].Data(),centNames[i].Data(),"Charged",sample),"READ");

    // if(i != 2 && i != 5)  file[i] = TFile::Open(Form("Result_%s_%s_%s_%s.root",mcType.Data(),corr[i].Data(),centNames[i].Data(),species.Data()),"READ");
    // else file[i] = TFile::Open(Form("Result_%s_%s_%s_%s.root",mcType.Data(),corr[i].Data(),centNames[i].Data(),"Charged"),"READ");
    // if(i == 2) file[i] = TFile::Open(Form("Result_%s_Nch_%s.root",corr[i].Data(),species.Data()),"READ");
    if(!file[i]) { printf("File %d not open \n",i); }
  }

 gStyle->SetOptStat(0);

 // Double_t etaMin = -1.6;
 // Double_t etaMax =  1.6;
 // Double_t etaMin[3] = {-5.0, 1.2,3.8};
 // Double_t etaMax[3] = {-1.1, 4., 7.};
 // Double_t etaMin[3] = {-5.0, 1.5, 4.6};
 // Double_t etaMax[3] = {-1.5, 3.8, 5.1};
 Double_t etaMin[3] = {-5.5, 1.4, 3.6};
 Double_t etaMax[3] = {-1.4, 4.0, 6.8};
 Int_t rebinDeta = 1;
 Int_t rebinDphi = 1;

 TH1D* hFinalV2[10] = {nullptr};
 for(Int_t p(0); p < 10; p++) hFinalV2[p] = new TH1D(Form("hFinalV2_%s_%d",species.Data(),p), "v_{2}; #it{p}_{T}; v_{2}", 13, binsV2);
 Double_t refFlow[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
 TString arrName[10] = {"TPC-FMDA", "TPC-FMDC", "FMDA-FMDC", "TPC-FMDA", "TPC-FMDC", "FMDA-FMDC", "TPC-FMDA", "TPC-FMDC", "FMDA-FMDC", ""};

 TH2D* hist[6];
 TCanvas* can[6];

 for(Int_t i(0); i < nt; i++){
   if(withAss && i == 1) continue;

   // if(i < 8 || i > 10) continue;
   // if(i > 3) continue;

   for(Int_t j(0); j < na; j++){
     // if(binst[i] < binsa[j]) continue;

     if(sameBin && (binst[i] != binsa[j])) continue;

     TLatex* latex = 0;
     latex = new TLatex();
     latex->SetTextSize(0.038);
     latex->SetTextFont(42);
     latex->SetTextAlign(21);
     latex->SetNDC();

     for(Int_t h(0); h < 6; h++){
       // if(h%3 != 0) continue; //testing
       if((h == 2 || h == 5) && i > 0) continue;
       hist[h] = (TH2D*)file[h]->Get(Form("dphi_%d_%d_%d",i,0,0));
       if(!hist[h]) { printf("Hist %d not open \n",h); }

       hist[h]->SetDirectory(0);
       hist[h]->GetXaxis()->SetTitleOffset(1.6);
       hist[h]->GetYaxis()->SetTitleOffset(1.6);
       hist[h]->GetZaxis()->SetTitleOffset(2.2);
       hist[h]->SetTitle(Form(";%s;%s;%s",kTitlePhi,kTitleEta,kCorrFuncTitle_Before));
       hist[h]->Rebin2D(rebinDphi,rebinDeta);
       hist[h]->Scale(1./rebinDphi/rebinDeta);

       if(i == 0){
         can[h] = new TCanvas(Form("cCan_%d",h),"cCent1",1000,1000);
         PadFor2DCorr();
         hist[h]->GetYaxis()->SetRangeUser(etaMin[h%3]+0.001,etaMax[h%3]-0.001);
         hist[h]->DrawCopy("surf1 FB");
         latex->DrawLatex(0.25,0.93,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
         latex->DrawLatex(0.25,0.86,centLegend[h].Data());
         TString textTr = Form("%.2f < p_{T}^{Tri} (GeV/c) < %.2f",binst[i],binst[i+1]);
         if(h != 2 && h != 5) latex->DrawLatex(0.78,0.86,textTr.Data());
         latex->DrawLatex(0.78,0.93,arrName[h].Data());
         // latex->DrawLatex(0.80,0.93,textTrg[h].Data());
         // latex->DrawLatex(0.80,0.86,textAss[h].Data());
         // gPad->Print(Form("%s/%s_%s_%s_%d.pdf",folder.Data(),species.Data(),centNames[h].Data(),corr[h].Data(),i));
       }


     }

     TH1D* proj=0x0;

     TH2D* histClone = nullptr;
     for(Int_t p(0); p < 6; p++){

       if((p == 2 || p == 5) && i > 0) continue;

       histClone = (TH2D*)hist[p]->Clone();
       if(!histClone) {printf("**** HistClone not found! **** \n"); return; }

       Int_t binEtaMin  = hist[p]->GetYaxis()->FindBin(etaMin[p%3] + 0.001);
       Int_t binEtaMax  = hist[p]->GetYaxis()->FindBin(etaMax[p%3] - 0.001);

       proj  = LinearFitProjectionX(histClone,Form("%s_proj1x", hist[p]->GetName()), binEtaMin, binEtaMax,0,0);

       proj->Scale(1./(binEtaMax - binEtaMin + 1));

       proj->SetStats(0);
       proj->GetYaxis()->SetNdivisions(505);
       proj->GetXaxis()->SetTitle("#Delta#varphi (rad)");
       proj->GetYaxis()->SetTitle(kProjYieldTitlePhi);
       proj->GetXaxis()->SetTitleOffset(1.1);
       proj->GetYaxis()->SetTitleOffset(1.1);
       proj->GetYaxis()->SetLabelSize(0.04);
       proj->GetXaxis()->SetLabelSize(0.04);
       proj->GetXaxis()->SetTitleSize(0.04);
       proj->GetYaxis()->SetTitleSize(0.04);
       proj->SetTitle("");
       proj->SetMarkerStyle(21);
       proj->SetMarkerSize(0.7);
       Float_t range = proj->GetMaximum()-proj->GetMinimum();
       proj->SetMaximum(proj->GetMinimum()+range*2);
       proj->SetMinimum(proj->GetMinimum()-range*0.1);

       if(p < 3){
         hminuit[p] = (TH1D*)proj->Clone(Form("hminuit_%d",p));
         continue;
       }

       counterCorr = p%3;

       hminuit_peri[counterCorr] = (TH1D*)proj->Clone(Form("hminuit_peri_%d",p));

       Double_t par[4],parerr[4];
       tempminuit23_1(par,parerr);

       hFinalV2[p]->SetBinContent(i+1, par[1]);
       hFinalV2[p]->SetBinError(i+1, parerr[1]);

       /*
       if(p == 1){
         periH = (TH1D*) proj->Clone();
         new TCanvas();
         proj->Draw("E0 X0");
         continue;
       }

       TF1 *f = new TF1("f",fun_template,-0.5*TMath::Pi()+1e-6,1.5*TMath::Pi()-1e-6,5);
       proj->Fit(f,"I0Q" ,"");

       TCanvas* c = new TCanvas("c2", "c2", 1200, 800);
       gPad->SetMargin(0.12,0.01,0.12,0.01);
       proj->GetYaxis()->SetTitleOffset(1.3);
       proj->DrawCopy("E0 X0");
       f->Draw("same");

       Double_t* pval = f->GetParameters();
       const Double_t* perr = f->GetParErrors();

       Float_t v0val = pval[0];
       Float_t v0err = TMath::Sqrt(perr[0]*perr[0]);

        Float_t chi2 = f->GetChisquare();
        Int_t ndf  = f->GetNDF();
        Printf("Chi2: %f ndf: %d chi2/ndf: %f",chi2,ndf,chi2/ndf);
        */

       //  bl->SetParameter(0,min);
       //  v1->SetParameters(pval[0],pval[1]);
       //  v2->SetParameters(pval[0],pval[2]);
       //  v3->SetParameters(pval[0],pval[3]);
       //
        // Float_t v1val = pval[1]/v0val;
        // Float_t v2val = pval[2]/v0val;
        // Float_t v3val = pval[3]/v0val;
        // Float_t v1err = v1val*sqrt(pow(v0err/v0val,2.)+pow(perr[1]/pval[1],2.));
        // Float_t v2err = v2val*sqrt(pow(v0err/v0val,2.)+pow(perr[2]/pval[2],2.));
        // Float_t v3err = v3val*sqrt(pow(v0err/v0val,2.)+pow(perr[3]/pval[3],2.));


        // printf("v2 val: %f \t v2 err: %f \n", v2val, v2err);
        // printf("sqrt v2 val: %f \t v2 err: %f \n", sqrt(v2val), v2err);

       //
       //  hFinalV2[p]->SetTitle(Form("v_{2}, %s",arrName[p].Data()));
       //
       //  if(sameBin) {
       //    hFinalV2[p]->SetBinContent(i+1, v2val);
       //    hFinalV2[p]->SetBinError(i+1, v2err);
       //
       //  }

        /*
        TCanvas* c = new TCanvas("c2", "c2", 1200, 800);
          gPad->SetMargin(0.12,0.01,0.12,0.01);
          proj->GetYaxis()->SetTitleOffset(1.3);
          proj->DrawCopy("E0 X0");
          vn->DrawCopy("same");
          bl->DrawCopy("same");
          v1->DrawCopy("same");
          v2->DrawCopy("same");
          v3->DrawCopy("same");

          TLegend* legend = new TLegend(0.48, 0.62, 0.98, 0.98);
          legend->SetFillColor(0);
          legend->SetBorderSize(0);
          legend->SetTextSize(0.035);

          legend->AddEntry(proj, Form("Data (%s)",species.Data()), "P");
          legend->AddEntry(vn, "a_{0} + #Sigma_{n=1}^{3} 2a_{n} cos(n#Delta#varphi)", "L");

          legend->AddEntry(v1, Form("a_{0} + 2a_{1} cos(#Delta#varphi), V_{1} #times10^{3}=%.2f#pm%.3f",v1val*1000,TMath::Abs(v1err*1000)), "L");
          legend->AddEntry(v2, Form("a_{0} + 2a_{2} cos(2#Delta#varphi), V_{2} #times10^{3}=%.2f#pm%.3f",v2val*1000,TMath::Abs(v2err*1000)), "L");
          legend->AddEntry(v3, Form("a_{0} + 2a_{3} cos(3#Delta#varphi), V_{3} #times10^{3}=%.2f#pm%.3f",v3val*1000,TMath::Abs(v3err*1000)), "L");

          legend->AddEntry(v1, "n=1", "L");
          legend->AddEntry(v2, "n=2", "L");
          legend->AddEntry(v3, "n=3", "L");


          legend->AddEntry(bl, "Baseline", "L");
          legend->Draw();

          latex->DrawLatex(0.3,0.93,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
          latex->DrawLatex(0.3,0.86,arrName[p].Data());
          latex->DrawLatex(0.3,0.79,textTrg.Data());
          latex->DrawLatex(0.3,0.72,textAss.Data());
          latex->DrawLatex(0.3,0.65,Form("%.1f < #Delta#eta < %.1f",etaMin,etaMax));

          gPad->Print(Form("%s/%s_fit_%d_%d_%d.pdf",folder.Data(),species.Data(),p,i,j));
          */

     } // end p


   } // end i
 } // end j

 Int_t harm = 2;
 TFile* fout = TFile::Open(Form("%s/v2_TF_%s_sample_%d.root",folder.Data(),species.Data(),sample),"RECREATE");

 for(Int_t p(0); p < 6; p++){
     // if(useAssRef){
     //   TFile* fInRef = TFile::Open(Form("plotsREF/v%d_TF_Charged.root",harm),"READ");
     //   TH1D* hInRef = (TH1D*)fInRef->Get(Form("hFinalV2_Charged_%d",p));
     //    if(!hInRef) { printf("Ref flow not loaded!\n"); return; }
     //    Double_t refContent = hInRef->GetBinContent(1);
     //    Double_t refError = hInRef->GetBinError(1);
     //    for(Int_t iBin(0); iBin < hFinalV2[p]->GetNbinsX()+2; iBin++){
     //      hFinalV2[p]->SetBinContent(iBin, hFinalV2[p]->GetBinContent(iBin)/refContent);
     //      hFinalV2[p]->SetBinError(iBin, sqrt( pow(hFinalV2[p]->GetBinError(iBin)/refContent,2.) +  pow(refError*hFinalV2[p]->GetBinContent(iBin),2.)/pow(refContent,4.) ) );
     //    }
     // }
     // else{
     //   for(Int_t iBin(0); iBin < hFinalV2[p]->GetNbinsX()+2; iBin++){
     //     if(hFinalV2[p]->GetBinContent(iBin) > 0.0) {
     //       hFinalV2[p]->SetBinContent(iBin, sqrt(hFinalV2[p]->GetBinContent(iBin)));
     //       hFinalV2[p]->SetBinError(iBin, sqrt(pow(hFinalV2[p]->GetBinError(iBin),2.)/hFinalV2[p]->GetBinContent(iBin))/2.);
     //     }
     //     else{
     //       hFinalV2[p]->SetBinContent(iBin, -10.0);
     //     }
     //   }
     // }
   // }
   TCanvas* cV2 = new TCanvas();
   hFinalV2[p]->Draw();
   fout->cd();
   hFinalV2[p]->Write();
   // cV2->SaveAs(Form("%s/v2_%s_%d_2.pdf",folder.Data(),species.Data(),p));
 }
 if(combine){
   hFinalV2[9]=(TH1D*)hFinalV2[0]->Clone("final");
   hFinalV2[9]->SetTitle("Combined, subtracted");
   Double_t vfmd = hFinalV2[5]->GetBinContent(1);
   Double_t vfmdEr = hFinalV2[5]->GetBinError(1);
   for(Int_t iBin(0); iBin < hFinalV2[0]->GetNbinsX()+2; iBin++){
     Double_t vfw = hFinalV2[3]->GetBinContent(iBin);
     Double_t vbw = hFinalV2[4]->GetBinContent(iBin);
     Double_t value = (vbw*vfw)/vfmd;

     Double_t vfwEr = hFinalV2[3]->GetBinError(iBin);
     Double_t vbwEr = hFinalV2[4]->GetBinError(iBin);

     if(value > 0.0) hFinalV2[9]->SetBinContent(iBin,sqrt(value));
     else hFinalV2[9]->SetBinContent(iBin,-9.9);

     Double_t error = (vfwEr*vfwEr*vbw)/(4.*vfw*vfmd) + (vbwEr*vbwEr*vfw)/(4.*vbw*vfmd) + (vfmdEr*vfmdEr*vfw*vbw)/(4.*vfmd*vfmd*vfmd);
     if(error > 0.0) hFinalV2[9]->SetBinError(iBin,sqrt(error));
     else hFinalV2[9]->SetBinError(iBin,9.9);
   }
   TCanvas* cV2 = new TCanvas();
   hFinalV2[9]->Draw();
   fout->cd();
   hFinalV2[9]->Write();
   // cV2->SaveAs(Form("%s/V%d_%s_%d.pdf",folder.Data(),harm,species.Data(),9));
   cV2->SaveAs(Form("%s/V%d_TF_%s_%d_sample_%d.pdf",folder.Data(),harm,species.Data(),9,sample));
 }

}


//==============================================================
void PadFor2DCorr(){
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.17);
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(55);
  gPad->SetPhi(45);
}
//==============================================================
void SymmetrizeTo0toPi(TH1D* h){
  for(Int_t k=1; k<=h->GetNbinsX(); k++){//symmetrize in 0-Pi
    Int_t ibin = -999.;
    Double_t phi = h->GetBinCenter(k);
    if      (phi<0)  ibin = h->FindBin(TMath::Abs(phi));
    else if (phi>pi) ibin = h->FindBin(2*pi-phi);
    if (ibin<0) continue;//here we are within 0 an Pi
    h->SetBinContent(ibin,(h->GetBinContent(ibin)+h->GetBinContent(k))/2);
    h->SetBinError(ibin,TMath::Sqrt(h->GetBinError(ibin)*h->GetBinError(ibin)+h->GetBinError(k)*h->GetBinError(k))/2);
    h->SetBinContent(k,0.);
    h->SetBinError(k,0.);
  }
}
//==============================================================
TH1D* LinearFitProjectionX(TH2D* h2, TString name, const Int_t bin1, const Int_t bin2, Bool_t UseStandardProj=kFALSE, Int_t debug=0){
  Float_t ymin = h2->GetYaxis()->GetBinLowEdge(bin1);
  Float_t ymax = h2->GetYaxis()->GetBinUpEdge(bin2);
  printf("Fit in %f %f\n",ymin,ymax);
  TH1D* h = h2->ProjectionX(name,bin1,bin2);
  TCanvas* cdebug = 0;
  TLatex* latex = 0;
  if (debug>1) {
   cdebug = new TCanvas("linproj","linproj",1900,1000);
   cdebug->Divide(h2->GetNbinsX()/3,3,0.001,0.001);
   latex = new TLatex();
   latex->SetTextSize(0.04);
   latex->SetTextFont(42);
   latex->SetTextAlign(11);
   latex->SetNDC();
  }
  Int_t binPi = h2->GetXaxis()->FindBin(0.5*pi - 0.001);

   for (Int_t ix=1;ix<=h2->GetNbinsX();ix++) {
    if(debug == -1 && ix < binPi){
      h->SetBinContent(ix,0.0);
      h->SetBinError(ix,0.0);
      continue;
    }

    TH1D* hbin = h2->ProjectionY(Form("hbin%i",ix),ix,ix);
    hbin->SetLineColor(kBlue);
    hbin->SetLineWidth(2);
    TF1* pol=0x0;
    if(gStudySystematic == k_parabolic_fit)pol = new TF1("pol","pol2",ymin,ymax);
    else if(gStudySystematic == k_constant_fit)pol = new TF1("pol","pol0",ymin,ymax);
    else pol = new TF1("pol","pol1",ymin,ymax);
    hbin->Fit(pol,"Q0","",ymin,ymax);
    Double_t value = pol->Integral(ymin,ymax)/hbin->GetXaxis()->GetBinWidth(1);
    Double_t error = pol->IntegralError(ymin,ymax)/hbin->GetXaxis()->GetBinWidth(1);
    // Standard Projection
    if (UseStandardProj) value = hbin->IntegralAndError(hbin->FindFixBin(ymin+0.001),hbin->FindFixBin(ymax-0.001),error);
    h->SetBinContent(ix,value);
    h->SetBinError(ix,error);

    if (cdebug) {
      cdebug->cd(ix);
      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetBottomMargin(0.11);
      hbin->SetMaximum(h2->GetMaximum()*1.1);
      hbin->SetMinimum(h2->GetMinimum()*0.98);
      hbin->DrawCopy();
      latex->DrawLatex(0.2,0.90,Form("#varphi bin = %i",ix));
      latex->DrawLatex(0.2,0.84,Form("#chi^{2}/ndf = %.2f/%i",pol->GetChisquare(),pol->GetNDF()));
      pol->SetLineWidth(1);
      pol->DrawCopy("same");
      value/=(bin2-bin1+1);
      error/=(bin2-bin1+1);
      latex->DrawLatex(0.2,0.78,Form("fit integral = %.4f #pm %.4f",value,error));
      value = hbin->IntegralAndError(hbin->FindFixBin(ymin+0.001),hbin->FindFixBin(ymax-0.001),error);
      value/=(bin2-bin1+1);
      error/=(bin2-bin1+1);
      latex->DrawLatex(0.2,0.72,Form("hist integral = %.4f #pm %.4f",value,error));
    }

    delete pol;
    delete hbin;
  }


  return h;
}
//==============================================================
Double_t fun_template(Double_t *x, Double_t *par)
{
  Int_t iBin = periH->GetXaxis()->FindBin(x[0]);
  Double_t periCont = periH->GetBinContent(iBin);
  return par[0]*periCont+par[1]*(1+2*par[2]*cos(1*x[0])+2*par[3]*cos(2*x[0])+2*par[4]*cos(3*x[0]));
}
//==============================================================
void minuitfcn23_1(int &npar, Double_t *gin, Double_t &ff, Double_t *par,int iflag){
  TH1D*h=(TH1D*)hminuit[counterCorr]->Clone("h");
  TH1D*hperi=(TH1D*)hminuit_peri[counterCorr]->Clone("hperi");
  Double_t f=par[0];
  Double_t gv=par[1];
  Double_t v2=par[2];
  Double_t v3=par[3];
  double lnQ = 0;
  for(int ibin=0; ibin<h->GetNbinsX(); ibin++){
    double x    = h->GetBinCenter(ibin+1);
    Double_t data=h->GetBinContent(ibin+1);
    Double_t exp=f*hperi->GetBinContent(ibin+1)+gv*(1.+2.*v2*cos(2.*x)+2.*v3*cos(3.*x));
    Double_t hoge=sqrt(hperi->GetBinError(ibin+1)*hperi->GetBinError(ibin+1)+h->GetBinError(ibin+1)*h->GetBinError(ibin+1));//sqrt(eybg[i]*eybg[i]+eysig[i]*eysig[i]);
    lnQ+=(data-exp)*(data-exp)/(hoge*hoge);
  }
  ff=lnQ;
}
//==============================================================
void tempminuit23_1(Double_t*fParamVal,Double_t*fParamErr)
{
  TFitter* minimizer = new TFitter(4);
  minimizer->SetFCN(minuitfcn23_1);
  minimizer->SetParameter(0,"Y", 0, 10, 0, 0);
  minimizer->SetParameter(1,"G", 0, 10, 0, 0);
  minimizer->SetParameter(2,"v2", 0, 1, 0, 0);
  minimizer->SetParameter(3,"v3", 0, 1, 0, 0);

  //  minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);

  //  double p0, ep0;
  double  p0 = minimizer->GetParameter(0);
  //double p1, ep1;
  double  p1 = minimizer->GetParameter(1);
  //  double p2, ep2;
  double p2 = minimizer->GetParameter(2);
  //  double p3, ep3;
  double  p3 = minimizer->GetParameter(3);

  fParamVal[0]=p0;//f
  fParamVal[1]=p2;//v2
  fParamVal[2]=p3;//v3
  fParamVal[3]=p1;//gv

  fParamErr[0]=minimizer->GetParError(0);
  fParamErr[1]=minimizer->GetParError(2);
  fParamErr[2]=minimizer->GetParError(3);
  fParamErr[3]=minimizer->GetParError(1);

  //  cout<<p0<<" "<<p1<<" "<<p2<<" "<<p3<<endl;
  cout << "Counter corr: " << counterCorr << "\n";
  cout << "F=" << fParamVal[0] << "+/-" << fParamErr[0] << "\n";
  cout << "v2=" << fParamVal[1] << "+/-" << fParamErr[1] << "\n";
  cout << "v3=" << fParamVal[2] << "+/-" << fParamErr[2] << "\n";
  cout << "gv=" << fParamVal[3] << "+/-" << fParamErr[3] << "\n";
}
//==============================================================
