void drawTest(TString species = "Charged", TString centrality = "cent"){
  // TFile *file = TFile::Open("AR_test.root","READ");
  // TFile *file = TFile::Open("/home/alidock/Downloads/AnalysisResults_DH.root","READ");
  TFile *file = TFile::Open("/home/alidock/output/pPb_LHC16q/train_3992/AnalysisResults.root","READ");
  // TFile *file = TFile::Open("/home/alidock/Downloads/AnalysisResults_2.root","READ");
  file->cd("CorrForFlow");
  TList* listIn = (TList*) gDirectory->Get(Form("Charged_Corr_%s",centrality.Data()));
  listIn->ls();

  AliTHn *contSE = nullptr;
  if(species == "Charged") contSE = (AliTHn*)listIn->FindObject("fhChargedSE");
  else if(species == "Pion" || species == "Kaon" || species == "Proton") contSE = (AliTHn*)listIn->FindObject(Form("fhPidSE_%s",species.Data()));
  else printf("***** This should not happen. ***** \n");
  contSE->FillParent();

  AliTHn *contME = nullptr;
  if(species == "Charged") contME = (AliTHn*)listIn->FindObject("fhChargedME");
  else if(species == "Pion" || species == "Kaon" || species == "Proton") contME = (AliTHn*)listIn->FindObject(Form("fhPidME_%s",species.Data()));
  else printf("***** This should not happen. ***** \n");
  contME->FillParent();

  AliCFGridSparse *sparSE = contSE->GetGrid(0);
  AliCFGridSparse *sparME = contME->GetGrid(0);

  TH2D* trig = nullptr;
  if(species == "Charged") trig = (TH2D*)listIn->FindObject("fhTrigTracks");
  else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH2D*)listIn->FindObject(Form("fhTrigTracks_%s",species.Data()));

  // TH1D* h1 = (TH1D*) sparSig->Project(1);
  // TH2D* h2 = (TH2D*) sparSE->Project(0,1);
  TH2D* h2 = (TH2D*) sparSE->Project(2,3);
  h2->Draw("colz");
  // h1->Draw("histo");

  TString outputC = "";
  if(centrality == "cent") outputC = "0_20";
  else outputC = "60_100";


  TFile *out = TFile::Open(Form("AR_%s_%s.root",outputC.Data(),species.Data()),"RECREATE");
  sparSE->Write();
  sparME->Write();
  trig->Write();


  return;
}
