void compToHEP(){
  // TFile *fileZ = TFile::Open("plots/v2.root","READ");
  TFile *fileZ = TFile::Open("plots/v2_Charged.root","READ");
  TFile *fileHEP = TFile::Open("HEPData-ins1242302-v1-root.root","READ");
  if(!fileZ || !fileHEP) {printf("Error with files \n"); return; }

  TDirectory* dir = fileHEP->GetDirectory("Table 1"); // 5
  // TDirectory* dir = fileHEP->GetDirectory("Table 17");
  TGraphAsymmErrors* gr = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
  gr->SetTitle("0-20%");
  // gr->SetTitle("60-100%");
  gr->GetXaxis()->SetRangeUser(0.2, 10.);

  TCanvas* c = new TCanvas();
  // gr->Draw();

  TH1D* h = (TH1D*)fileZ->Get("hFinalV2_Charged_2");
  TGraphErrors* grH = new TGraphErrors(h);

  grH->SetLineColor(kRed);
  grH->SetTitle("(0-20%)-(60-100%)");
  grH->GetYaxis()->SetRangeUser(0.0, 0.2);
  grH->GetXaxis()->SetRangeUser(0.2, 4.0);

  // grH->Draw("same");
  grH->Draw();
  gr->Draw("same");

  c->SaveAs("new_sub.pdf");


  return;
}
