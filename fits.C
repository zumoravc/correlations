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

Double_t binst[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
// Double_t binst[]={1.0, 2., 1.0, 2.0, 4.0};
const Int_t nt = sizeof(binst) / sizeof(Double_t) - 1;

const Double_t binsa[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
const Int_t na = sizeof(binsa) / sizeof(Double_t) - 1;

// Double_t binsV2[9] = {0.5, 0.7, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0};
// Double_t binsV2[7] = {1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0};
Double_t binsV2[22] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};

Bool_t withAss = kFALSE;
Bool_t sameBin = kTRUE;

void fits(        TString species = "Pion",
                  TString sStyle = "CC",
                  TString sEst = "V0A",
                  Int_t c1 = 0,
                  Int_t c2 = 1,
                  Double_t* res = NULL)
{


 TFile *file1 = TFile::Open(Form("Result_TPC_0_20_%s.root",species.Data()),"READ");
 TFile *file2 = TFile::Open(Form("Result_TPC_60_100_%s.root",species.Data()),"READ");
 if(!file1 || !file2) { printf("%s\n"); }

 gStyle->SetOptStat(0);

 Double_t etaMin = -1.6;
 Double_t etaMax =  1.6;
 Int_t rebinDeta = 1;
 Int_t rebinDphi = 1;
 Float_t exclusionMin = -0.8;
 Float_t exclusionMax = 0.8;

 TH1D* hFinalV2[3] = {nullptr};
 for(Int_t p(0); p < 3; p++) hFinalV2[p] = new TH1D(Form("hFinalV2_%s_%d",species.Data(),p), "v_{2}; #it{p}_{T}; v_{2}", 21, binsV2);
 Double_t refFlow[3] = {0.0, 0.0, 0.0};

 for(Int_t i(0); i < nt; i++){
   if(withAss && i == 1) continue;

   for(Int_t j(0); j < na; j++){
     if(binst[i] < binsa[j]) continue;

     if(sameBin && (binst[i] != binsa[j])) continue;

     TH2D *hist1 = (TH2D*)file1->Get(Form("dphi_%d_%d_%d",i,j,c1));
     hist1->SetDirectory(0);

     TH2D *hist2 = (TH2D*)file2->Get(Form("dphi_%d_%d_%d",i,j,c2));
     hist2->SetDirectory(0);

     hist1->GetXaxis()->SetTitleOffset(1.6);
     hist1->GetYaxis()->SetTitleOffset(1.6);
     hist1->GetZaxis()->SetTitleOffset(2.2);
     hist2->GetXaxis()->SetTitleOffset(1.6);
     hist2->GetYaxis()->SetTitleOffset(1.6);
     hist2->GetZaxis()->SetTitleOffset(2.2);

     TString label1(hist1->GetTitle());
     TString label2(hist2->GetTitle());
     TObjArray* objArray1 = label1.Tokenize(",~");
     TObjArray* objArray2 = label2.Tokenize(",~");
     printf("label1=%s\n",label1.Data());
     printf("label2=%s\n",label2.Data());

     TString textTrg = objArray1->At(0)->GetName();
     TString textAss = objArray1->At(1)->GetName();
     TString cent11 = objArray1->At(2)->GetName();
     TString cent12 = objArray1->At(3)->GetName();
     TString cent21 = objArray2->At(2)->GetName();
     TString cent22 = objArray2->At(3)->GetName();
     cent11.ReplaceAll(" ","");
     cent21.ReplaceAll(" ","");
     TString textCent1   = Form("%s-%s", cent11.Data(),cent12.Data());
     TString textCent2   = Form("%s-%s", cent21.Data(),cent22.Data());
     TString textCentDif = Form("(%s) - (%s)", textCent1.Data(),textCent2.Data());
     hist1->SetTitle(Form(";%s;%s;%s",kTitlePhi,kTitleEta,kCorrFuncTitle_Before));
     hist2->SetTitle(Form(";%s;%s;%s",kTitlePhi,kTitleEta,kCorrFuncTitle_Before));

     TString arrName[3] = {textCent1, textCent2, textCentDif};

     hist1->Rebin2D(rebinDphi,rebinDeta);
     hist2->Rebin2D(rebinDphi,rebinDeta);
     hist1->Scale(1./rebinDphi/rebinDeta);
     hist2->Scale(1./rebinDphi/rebinDeta);


     TLatex* latex = 0;
     latex = new TLatex();
     latex->SetTextSize(0.04);
     latex->SetTextFont(42);
     latex->SetTextAlign(21);
     latex->SetNDC();


      TCanvas* cCent1 = new TCanvas("cCent1","cCent1",1000,1000);
      PadFor2DCorr();
      hist1->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
    //    hist1->GetZaxis()->SetRangeUser(3.15,3.4);
        //    hist1->SetMaximum(1.5*hist1->GetBinContent(hist1->FindBin(3.14159,0))-0.5*hist1->GetMinimum());
        hist1->DrawCopy("surf1 FB");
        latex->DrawLatex(0.25,0.93,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
        latex->DrawLatex(0.25,0.86,textCent1.Data());
        latex->DrawLatex(0.80,0.93,textTrg.Data());
        latex->DrawLatex(0.80,0.86,textAss.Data());
        // gPad->Print("0_20.pdf");
        gPad->Print(Form("plots/%s_0_20_%d_%d.pdf",species.Data(),i,j));
        //gPad->Print(Form("cent1%i%s",i,".pdf"));

        TCanvas* cCent2 = new TCanvas("cCent2","cCent2",1000,1000);
        PadFor2DCorr();
        hist2->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
        //    hist2->SetMaximum(1.5*hist2->GetBinContent(hist2->FindBin(3.14159,0))-0.5*hist2->GetMinimum());
        hist2->DrawCopy("surf1 FB");
        latex->DrawLatex(0.25,0.93,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
        latex->DrawLatex(0.25,0.86,textCent2.Data());
        latex->DrawLatex(0.80,0.93,textTrg.Data());
        latex->DrawLatex(0.80,0.86,textAss.Data());
        // gPad->Print("60_100.pdf");
        gPad->Print(Form("plots/%s_60_100_%d_%d.pdf",species.Data(),i,j));
        //gPad->Print(Form("cent2%i%s",i,".pdf"));

     Float_t phiBaselineMin_H = pi/2 - 0.2 + 0.0001;
     Float_t phiBaselineMin = 0. - 0.01 + 0.0001;
     Float_t phiBaselineMax_H = pi/2 + 0.2 - 0.0001;
     Float_t phiBaselineMax = 0. + 0.01 - 0.0001;
     Int_t binEtaMin  = hist1->GetYaxis()->FindBin(etaMin + 0.001);
     Int_t binEtaMax  = hist1->GetYaxis()->FindBin(etaMax - 0.001);

     Int_t binPhiBaselineMin = hist1->GetXaxis()->FindBin(phiBaselineMin);
     Int_t binPhiBaselineMax = hist1->GetXaxis()->FindBin(phiBaselineMax);

     Int_t binPhiBaselineMin_H = hist1->GetXaxis()->FindBin(phiBaselineMin_H);
     Int_t binPhiBaselineMax_H = hist1->GetXaxis()->FindBin(phiBaselineMax_H);

     Double_t baseLineE;
     Float_t baseLine = hist1->IntegralAndError(binPhiBaselineMin_H, binPhiBaselineMax_H, binEtaMin, binEtaMax, baseLineE);
     baseLine  /= (binPhiBaselineMax_H - binPhiBaselineMin_H + 1) * (binEtaMax - binEtaMin + 1);
     baseLineE /= (binPhiBaselineMax_H - binPhiBaselineMin_H + 1) * (binEtaMax - binEtaMin + 1);

     Double_t baseLineE2;
     Float_t baseLine2 = hist2->IntegralAndError(binPhiBaselineMin, binPhiBaselineMax, binEtaMin, binEtaMax, baseLineE2);
     baseLine2  /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);
     baseLineE2 /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);

     // hist1->Add(hist2, -1);

     TCanvas* cSubt = new TCanvas("cSubt", "cSubt",1000,1000);
     PadFor2DCorr();
     TH2D* histSub = (TH2D*) hist1->Clone("histSub");
     histSub->Add(hist2, -1);
     histSub->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
     histSub->DrawCopy("surf1 FB");
     latex->DrawLatex(0.25,0.93,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
     latex->DrawLatex(0.25,0.86,textCentDif.Data());
     latex->DrawLatex(0.80,0.93,textTrg.Data());
     latex->DrawLatex(0.80,0.86,textAss.Data());
     // gPad->Print("subt.pdf");
     gPad->Print(Form("plots/%s_subt_%d_%d.pdf",species.Data(),i,j));

     TH1D* proj=0x0;
     TH1D* proj2=0x0;

     Int_t binExclusionMin  = hist1->GetYaxis()->FindBin(exclusionMin - 0.001);
     Int_t binExclusionMax  = hist1->GetYaxis()->FindBin(exclusionMax + 0.001);

     Double_t minPer = 0.0;
     for(Int_t p(0); p < 3; p++){
       if(p == 0){
         proj  = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,0,0);
         proj2 = LinearFitProjectionX(hist1,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,0,0);
       }
       else if(p == 1){
         proj  = LinearFitProjectionX(hist2,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,0,0);
         proj2 = LinearFitProjectionX(hist2,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,0,0);
       }
       else{
         proj  = LinearFitProjectionX(histSub,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,0,0);
         proj2 = LinearFitProjectionX(histSub,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,0,0);
       }

       proj->Add(proj2, 1);
       proj->Scale(1./(binEtaMax - binExclusionMax + 1 + binExclusionMin - binEtaMin + 1));


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

       TF1 *vn = new TF1("vn","[0]+2*[1]*cos(1*x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)",-5,5);

       TF1* v1 = new TF1("v1","[0]+2*[1]*cos(1*x)",-5,5);
       TF1* v2 = new TF1("v2","[0]+2*[1]*cos(2*x)",-5,5);
       TF1* v3 = new TF1("v3","[0]+2*[1]*cos(3*x)",-5,5);
       TF1* bl = new TF1("bl","[0]",-5,5);
       vn->SetLineStyle(kSolid);       vn->SetLineColor(kBlack);
       bl->SetLineStyle(kSolid);       bl->SetLineColor(kBlue);
       v1->SetLineStyle(kDotted);      v1->SetLineColor(kBlue);
       v2->SetLineStyle(kDashed);      v2->SetLineColor(kRed);
       v3->SetLineStyle(kDashDotted);  v3->SetLineColor(kMagenta);

       if (gStudySystematic==k_no_v3) vn->FixParameter(3,0);
       proj->Fit(vn,"I0Q" ,"");

       Double_t* pval = vn->GetParameters();
       const Double_t* perr = vn->GetParErrors();
       Float_t min = vn->GetMinimum();

       if(p < 2) {
         baseLine2 = 0.0;
         baseLineE2 = 0.0;
       }

       Float_t v0val = baseLine2+pval[0];
       Float_t v0err = TMath::Sqrt(baseLineE2*baseLineE2+perr[0]*perr[0]);

        Float_t chi2 = vn->GetChisquare();
        Int_t ndf  = vn->GetNDF();
        Printf("Chi2: %f ndf: %d chi2/ndf: %f Min: %f",chi2,ndf,chi2/ndf, min);

        bl->SetParameter(0,min);
        v1->SetParameters(pval[0],pval[1]);
        v2->SetParameters(pval[0],pval[2]);
        v3->SetParameters(pval[0],pval[3]);

        Float_t v1val = pval[1]/v0val;
        Float_t v2val = pval[2]/v0val;
        Float_t v3val = pval[3]/v0val;
        Float_t v1err = v1val*sqrt(pow(v0err/v0val,2.)+pow(perr[1]/pval[1],2.));
        Float_t v2err = v2val*sqrt(pow(v0err/v0val,2.)+pow(perr[2]/pval[2],2.));
        Float_t v3err = v3val*sqrt(pow(v0err/v0val,2.)+pow(perr[3]/pval[3],2.));

        Printf("Baseline: %f +- %f; v1 = %f +- %f; v2 = %f +- %f; v3 = %f +- %f \n", baseLine, baseLineE, v1val, v1err, v2val, v2err, v3val, v3err);


        if(!sameBin && i == 0 && j == 0) {
          // if(v2val > 0.0) refFlow[p] = sqrt(pval[2]/(pval[0]+minPer));
          if(v2val > 0.0) refFlow[p] = sqrt(pval[2]/(min));
          else printf(" ******** Problem with ref. flow value ******** v2: %f \n", v2val);
        }

        if(withAss && i > 1){
          hFinalV2[p]->SetBinContent(i-1, (pval[2]/(min))/refFlow[p]);
        }

        hFinalV2[p]->SetTitle(Form("v_{2}, %s",arrName[p].Data()));

        if(sameBin) {
          if(v2val > 0.0){
            hFinalV2[p]->SetBinContent(i+1, sqrt(v2val));
            hFinalV2[p]->SetBinError(i+1, v2err/(2.0*sqrt(v2val)));
          }
          else{
            hFinalV2[p]->SetBinContent(i+1, -10.0);
            hFinalV2[p]->SetBinError(i+1, 1.0);
          }
        }

        if(p==1) minPer = pval[0];
        Printf("Baseline2: %f +- %f; a0 from peri: %f  \n", baseLine2, baseLineE2, minPer);

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

          gPad->Print(Form("plots/%s_fit_%d_%d_%d.pdf",species.Data(),p,i,j));

          if (res) {
           res[0] = baseLine;
           res[1] = v1val;
           res[2] = v2val;
           res[3] = v3val;
           res[4] = baseLineE;
           res[5] = v1err;
           res[6] = v2err;
           res[7] = v3err;
           res[8] = chi2;
           res[9] = ndf;
         }

     } // end p

     // if(gStudySystematic == k_StandardProj)
     // {
     //  proj  = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,1,0);
     //  proj2 = LinearFitProjectionX(hist1,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,1,0);
     // }
     // else
     // {
     //  proj  = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,0,0);
     //  proj2 = LinearFitProjectionX(hist1,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,0,0);
     // }




   } // end i
 } // end j

 TFile* fout = TFile::Open(Form("plots/v2_%s.root",species.Data()),"RECREATE");

 TFile* fInRef = TFile::Open("plots/v2_Charged.root","READ");

 for(Int_t p(0); p < 3; p++){
   if(species != "Charged"){
     TH1D* hInRef = (TH1D*)fInRef->Get(Form("hFinalV2_Charged_%d",p));
     if(!hInRef) { printf("Ref flow not loaded!\n"); return; }
     if(hInRef->GetNbinsX() != hFinalV2[p]->GetNbinsX()) { printf("Dif. N of bins!\n"); return; }
     for(Int_t iBin(0); iBin < hFinalV2[p]->GetNbinsX()+2; iBin++){
       if(hFinalV2[p]->GetBinContent(iBin) > 0.0) {
         hInRef->SetBinContent(iBin, hInRef->GetBinContent(iBin)/hFinalV2[p]->GetBinContent(iBin));
         // hInRef->SetBinError(iBin, hInRef->GetBinContent(iBin)/hFinalV2[p]->GetBinContent(iBin));
       }
       else{
         hInRef->SetBinContent(iBin, -10.0);
       }
     }
   }
   TCanvas* cV2 = new TCanvas();
   hFinalV2[p]->Draw();
   hFinalV2[p]->Write();
   cV2->SaveAs(Form("plots/v2_%s_%d.pdf",species.Data(),p));
 }

}



void PadFor2DCorr(){
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.17);
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(55);
  gPad->SetPhi(45);
}
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
   for (Int_t ix=1;ix<=h2->GetNbinsX();ix++) {
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
