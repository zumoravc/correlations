#ifdef __CLING__
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCorrForFlow.cxx"
#endif
#include "AliMultSelectionTask.h"


void runTask()
{
    Bool_t local = 1;
    Bool_t gridTest = 0;

    #if !defined (__CINT__) || defined (__CLING__)
      gInterpreter->ProcessLine(".include $ROOTSYS/include");
      gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
      gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    #else
      gROOT->ProcessLine(".include $ROOTSYS/include");
      gROOT->ProcessLine(".include $ALICE_ROOT/include");
      gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    #endif

    AliAnalysisManager *mgr = new AliAnalysisManager("AliAnalysisTaskCorrForFlow");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    AliMultSelectionTask* taskMultSelection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);
    // taskMultSelection->SetSelectedTriggerClass(AliVEvent::kHighMultV0);

    AliAnalysisTaskPIDResponse* taskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"))));

#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskCorrForFlow.cxx+g");
    AliAnalysisTaskCorrForFlow *task = reinterpret_cast<AliAnalysisTaskCorrForFlow*>(gInterpreter->ExecuteMacro("AddCorrForFlowTask.C(\"SomeName\", \"suffixTestik\")"));
#else
    gROOT->LoadMacro("AliAnalysisTaskCorrForFlow.cxx++g");
    gROOT->LoadMacro("AddCorrForFlowTask.C");
    AliAnalysisTaskCorrForFlow *task = AddCorrForFlowTask();
#endif

    // task->SetTrigger(AliVEvent::kHighMultV0);
    // task->SetTrigger(AliVEvent::kINT7);
    // task->AddFlowTask({2,-2});
    // // task->AddFlowTask({2,-2});
    // task->AddFlowTask({2,2,-2,-2});
    // task->SetMaxPowersVector({5,0,4,0,3});
    // task->SetNofSamples(5);
    task->SetCentrality("V0A", 0., 20.);
    task->SetPtBins({0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0});
    task->SetPtBinsAss({0.5, 1.0, 1.5, 2.0, 3.0});
    task->SetPtRangeTrig(0.2, 10.0);
    task->SetPtRangeAss(0.5, 3.0);
    task->SetDoPID();

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        TChain* chain = new TChain("aodTree");
        chain->Add("AliAOD.root");
        mgr->StartAnalysis("local", chain);
    } else {
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // alienHandler->SetAdditionalLibs("AliAnalysisTaskCorrForFlow.h AliAnalysisTaskCorrForFlow.cxx");
        alienHandler->SetAdditionalLibs("PWGCFFLOWGF.par libPWGEMCALbase.so");
        alienHandler->SetAnalysisSource("AliAnalysisTaskCorrForFlow.cxx");
        alienHandler->SetAliPhysicsVersion("vAN-20210628_ROOT6-1");
        // alienHandler->SetAPIVersion("V1.1x");
        alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");
        alienHandler->SetDataPattern("pass1_CENT_woSDD/AOD190/*/AliAOD.root");
        alienHandler->SetRunPrefix("000");
        alienHandler->AddRunNumber(265338);
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("Corr13.sh");
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("Corr13.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        alienHandler->SetGridWorkingDir("par13");
        alienHandler->SetGridOutputDir("output");

        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            alienHandler->SetNtestFiles(1);
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            alienHandler->SetRunMode("full");
            // alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}
