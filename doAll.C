{

    gROOT->ProcessLine(".L CORE/CMS3_CORE.so");
    gROOT->ProcessLine(".L ScanChain.C+");

    
    TChain *ch_tt = new TChain("Events","tt");
    ch_tt->Add("/hadoop/cms/store/group/snt/run2_mc2017/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05/merged_ntuple_1.root");
    ScanChain(ch_tt, "tt", true);

}

