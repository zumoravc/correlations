void loop(TString corr, TString species, TString centralityOutput);

void check(){
  // loop("", "K0s", "60_100");
  // loop("", "Lambda", "60_100");

  TString species = "Charged";
  TString cent = "HM";
  // TString cent = "0_20";
  // TString flag = "PVz";

  // loop("gen", species, cent);
  // loop("reco", species, cent);

  loop("", species, cent);
  // loop(flag, species, cent);
  // loop("2", species, cent);
  // loop("3", species, cent);
  // loop("FMDAFMDC", species, cent);
  // loop("TPCFMDA", species, cent);
  // loop("TPCFMDC", species, cent);



  TString specList[3] = {"Pion", "Kaon", "Proton"};
  // TString specList[2] = {"K0s", "Lambda"};
  for(Int_t i(0); i < 3; i++){
    // loop("gen", specList[i], cent);
    // loop("reco", specList[i], cent);
    // loop("TPCFMDC", specList[i], cent);
    // loop("TPCFMDC", specList[i], cent);
    loop("", specList[i], cent);
    // loop(flag, specList[i], cent);
  }


}



void loop(TString corr, TString species, TString centralityOutput){
  TString path = "/Users/zuzana/Mirror/output/pp_LHC16/train_6390";
  // TString path = "/home/alidock/output/pp_LHC18/FMD_train_628_56/HM";
  // TString path = "/Users/zuzana/Mirror/output/pPb_LHC16q/train_4228";
  // TString path = "/Users/zuzana/Mirror/output/pPb_MC/train_1330";
  // TString path = "/Users/zuzana/Downloads";
  // TString path = "/home/alidock/correlations/taskFMD";
  TFile *file = TFile::Open(Form("%s/AnalysisResults.root",path.Data()),"READ");
  // TFile *file = TFile::Open(Form("%s/grid.root",path.Data()),"READ");
  // TFile *file = TFile::Open(Form("%s/AR.root",path.Data()),"READ");
  file->cd("CorrForFlow");
  file->ls();
  TList* listIn = (TList*) gDirectory->Get(Form("Charged_Corr_%s",corr.Data()));
  // TList* listIn = (TList*) gDirectory->Get(Form("Charged_%s_",corr.Data()));
  // TList* listIn = (TList*) gDirectory->Get(Form("Charged_%s_%s",corr.Data(),centralityOutput.Data()));
  listIn->ls();

  //testing + checking
  /*
  AliTHn *contSE = (AliTHn*)listIn->FindObject("fhChargedSE");
  // AliTHn *contSE = (AliTHn*)listIn->FindObject("fhPidSE_Proton");
  if(!contSE) {printf("Problem with contSE \n"); return; }
  contSE->FillParent();
  AliCFGridSparse *sparSE = contSE->GetGrid(0);
  TH2D* h2 = (TH2D*) sparSE->Project(0,1);
  new TCanvas();
  h2->Draw("surf1");
  TH2D* h22 = (TH2D*) sparSE->Project(0,2);
  new TCanvas();
  h22->Draw("surf1");

  AliTHn *contME = (AliTHn*)listIn->FindObject("fhChargedME");
  if(!contME) {printf("Problem with contSE \n"); return; }
  contME->FillParent();
  AliCFGridSparse *sparME = contME->GetGrid(0);
  TH2D* h2M = (TH2D*) sparME->Project(0,1);
  new TCanvas();
  h2M->Draw("surf1");
  TH2D* h22M = (TH2D*) sparME->Project(0,2);
  new TCanvas();
  h22M->Draw("surf1");

  TH2D* trig = (TH2D*)listIn->FindObject("fhTrigTracks");
  if(!trig) {printf("Problem with trig \n"); return; }
  new TCanvas();
  trig->Draw("colz");
  */


  AliTHn *contSE = nullptr;
  if(species == "Charged") contSE = (AliTHn*)listIn->FindObject("fhChargedSE");
  else if(species == "Pion" || species == "Kaon" || species == "Proton" || species == "K0s" || species == "Lambda") contSE = (AliTHn*)listIn->FindObject(Form("fhPidSE_%s",species.Data()));
  else printf("***** This should not happen. ***** \n");
  contSE->FillParent();

  AliTHn *contME = nullptr;
  if(species == "Charged") contME = (AliTHn*)listIn->FindObject("fhChargedME");
  else if(species == "Pion" || species == "Kaon" || species == "Proton" || species == "K0s" || species == "Lambda") contME = (AliTHn*)listIn->FindObject(Form("fhPidME_%s",species.Data()));
  else printf("***** This should not happen. ***** \n");
  contME->FillParent();

  AliCFGridSparse *sparSE = contSE->GetGrid(0);
  AliCFGridSparse *sparME = contME->GetGrid(0);

  TH3D* trig = nullptr;
  // if(species == "Charged") trig = (TH3D*)listIn->FindObject("fhTrigTracks");
  // else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH3D*)listIn->FindObject(Form("fhTrigTracks_%s",species.Data()));

  AliTHn* contTrig = nullptr;
  AliCFGridSparse *sparsTrig = nullptr;
  // if(species == "K0s" || species == "Lambda") {
  if(kTRUE) {
    if(species == "Charged") contTrig = (AliTHn*)listIn->FindObject(Form("fhTrigTracks"));
    else contTrig = (AliTHn*)listIn->FindObject(Form("fhTrigTracks_%s",species.Data()));
    contTrig->FillParent();
    sparsTrig = contTrig->GetGrid(0);
  }

  // TH2D* trig = nullptr;
  // if(species == "Charged") trig = (TH2D*)listIn->FindObject("fhTrigTracks");
  // else if(species == "Pion" || species == "Kaon" || species == "Proton") trig = (TH2D*)listIn->FindObject(Form("fhTrigTracks_%s",species.Data()));

  // if(centralityOutput == "") centralityOutput = centrality;

  TFile *out = TFile::Open(Form("%s/AR_%s_%s.root",path.Data(),centralityOutput.Data(),species.Data()),"RECREATE");
  // TFile *out = TFile::Open(Form("%s/AR_%s_%s_%s.root",path.Data(),corr.Data(),centralityOutput.Data(),species.Data()),"RECREATE");
  sparSE->Write();
  sparME->Write();
  if(sparsTrig) sparsTrig->Write();
  if(trig) trig->Write();


  return;
}
