#pragma GCC diagnostic ignored "-Wsign-compare"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1.h"
#include "TChain.h"

#include <math.h>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

#include "CORE/Tools/JetCorrector.h"
#include "CORE/CMS3.h"
#include "CORE/JetSelections.h"
#include "CORE/IsolationTools.cc"
#include "CORE/LeptonSelections.cc"
#include "CORE/ElectronSelections.cc"
#include "CORE/SSSelections.cc"
#include "CORE/VertexSelections.cc"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"
#include "CORE/Config.h"

#include "tqdm.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;

bool STOP_REQUESTED = false;

int is_truth_matched(int id, int idx) {
    // charge flips fall into -1
    Lep lep = Lep(id,idx);
    int mid = lepMotherID_v2(lep).first;
    if (mid == 1) return 1;
    else if (mid <= 0) return 0;
    else return -1;
}

int get_mother_id(int id, int idx) {
    Lep lep = Lep(id,idx);
    return lepMotherID_v2(lep).first;
}

int ScanChain(TChain *ch, TString suffix, bool isMC=true){

    signal(SIGINT, [](int){
            cout << "SIGINT Caught, stopping after current event" << endl;
            STOP_REQUESTED=true;
            });
    STOP_REQUESTED=false;

    // TH1F * h_met = new TH1F("met", "met", 50, 0, 300);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);

    tqdm bar;

    DatasetInfoFromFile df;
    df.loadFromFile("CORE/Tools/datasetinfo/scale1fbs.txt");

    createAndInitMVA("CORE", true, true, 80);

    gconf.year = 2017;
    gconf.SS_innerlayers = 0;
    gconf.ea_version = 4;
    gconf.btag_disc_wp = 0.4941;
    gconf.multiiso_el_minireliso = 0.09;
    gconf.multiiso_el_ptratio = 0.85;
    gconf.multiiso_el_ptrel = 9.2;
    gconf.multiiso_mu_minireliso = 0.12;
    gconf.multiiso_mu_ptratio = 0.80;
    gconf.multiiso_mu_ptrel = 7.5;

    //JEC files -- 25 ns MC
    std::string jecEraMC = "Fall17_17Nov2017_V6";
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL1L2L3;
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L1FastJet_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L3Absolute_AK4PFchs.txt");
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL2L3;
    jetcorr_filenames_25ns_MC_pfL2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L3Absolute_AK4PFchs.txt");
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL1;
    jetcorr_filenames_25ns_MC_pfL1.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L1FastJet_AK4PFchs.txt");
    gconf.jet_corrector_L1 = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL1);
    gconf.jet_corrector_L2L3 = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL2L3);

    float lumiAG = 41.3;

    TFile *out_file = new TFile("output_"+suffix+".root", "RECREATE");
    out_file->cd();
    TTree* out_tree = new TTree("t", "tuneiso");

    float eta;
    float phi;
    float ht;
    float miniiso;
    float pt;
    float ptratio;
    float ptrel;
    int nvtx;
    int id;
    int sample;
    int njets;
    int ssid;
    int truth;
    int mother_id;
    float pucorr;
    float weight;
    float rho;
    float miniiso_ch;
    float miniiso_nh;
    float miniiso_em;
    float miniiso_raw;
    float miniiso_ncorr;
    float disc;
    float nmiss;
    float hovere;
    float met;
    int ngenlepw;
    int pass_iso;
    int pass_id;

    out_tree->Branch("eta", &eta);
    out_tree->Branch("phi", &phi);
    out_tree->Branch("ht", &ht);
    out_tree->Branch("miniiso", &miniiso);
    out_tree->Branch("pt", &pt);
    out_tree->Branch("ptratio", &ptratio);
    out_tree->Branch("ptrel", &ptrel);
    out_tree->Branch("id", &id);
    out_tree->Branch("sample", &sample);
    out_tree->Branch("njets", &njets);
    out_tree->Branch("ssid", &ssid);
    out_tree->Branch("truth", &truth);
    out_tree->Branch("mother_id", &mother_id);
    out_tree->Branch("nvtx", &nvtx);
    out_tree->Branch("weight", &weight);
    out_tree->Branch("pucorr", &pucorr);
    out_tree->Branch("rho", &rho);
    out_tree->Branch("miniiso_ch", &miniiso_ch);
    out_tree->Branch("miniiso_nh", &miniiso_nh);
    out_tree->Branch("miniiso_em", &miniiso_em);
    out_tree->Branch("miniiso_raw", &miniiso_raw);
    out_tree->Branch("miniiso_ncorr", &miniiso_ncorr);
    out_tree->Branch("ngenlepw", &ngenlepw);
    out_tree->Branch("disc", &disc);
    out_tree->Branch("met", &met);
    out_tree->Branch("nmiss", &nmiss);
    out_tree->Branch("pass_iso", &pass_iso);
    out_tree->Branch("pass_id", &pass_id);
    out_tree->Branch("hovere", &hovere);

    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        if (STOP_REQUESTED) break;

        TString name = ch->GetTitle();
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        cms3.Init(tree);

        TString filename(currentFile->GetTitle());

        auto tokens = filename.Tokenize("/");
        auto basename = ((TObjString*)(tokens->At(tokens->GetEntries()-1)))->String().Data();
        bar.set_label(basename);

        sample = 0;
        if (filename.Contains("TTTT")) sample = 1;

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
            if (STOP_REQUESTED) break;

            cms3.GetEntry(event);
            nEventsTotal++;

            // if (nEventsTotal > 200000) break;

            // CMS3::progress(nEventsTotal, nEventsChain);
            bar.progress(nEventsTotal, nEventsChain);

            // rho = evt_fixgridfastjet_centralneutral_rho();
            rho = evt_fixgridfastjet_all_rho();
            met = evt_pfmet();

            // float sgnMCweight = ((tas::genps_weight() > 0) - (tas::genps_weight() < 0));
            // weight = lumiAG*sgnMCweight*df.getScale1fbFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
            // weight *= getTruePUw(tas::puInfo_trueNumInteractions()[0]);
            weight = 1.;

            for (int imu = 0; imu < mus_p4().size(); imu++) {
                pt = mus_p4()[imu].pt();
                if (pt < 15.) continue;

                eta = mus_p4()[imu].eta();
                if (fabs(eta) > 2.4) continue;
                id = -13*mus_charge()[imu];

                phi = mus_p4()[imu].phi();

                // ssid = isGoodLeptonNoIso(id,imu);
                // if (!ssid) continue;

                if (isMC) {
                    truth = is_truth_matched(id,imu);
                    mother_id = get_mother_id(id,imu);
                }
                miniiso = muMiniRelIsoCMS3_EA(imu, gconf.ea_version);
                ptrel = getPtRel(id, imu, true, ssWhichCorr);
                ptratio = pt/closestJet(mus_p4()[imu], 0.4, 3.0, ssWhichCorr).pt();
                pucorr = rho * muEA03(imu, gconf.ea_version) * pow(getMiniDR(pt)/0.3,2);

                miniiso_ch = mus_miniIso_ch()[imu] / pt;
                miniiso_nh = mus_miniIso_nh()[imu] / pt;
                miniiso_em = mus_miniIso_em()[imu] / pt;
                miniiso_raw = miniiso_ch + miniiso_nh + miniiso_em;
                miniiso_ncorr = miniiso_nh + miniiso_em - pucorr/pt;

                // pass_iso = (miniiso<0.12) && (ptratio>0.80 || ptrel>7.5);
                pass_id = isGoodLeptonNoIso(id,imu);
                pass_iso = isGoodLepton(id,imu);

                hovere = -1.;

                nmiss = 0;

                out_tree->Fill();

            }

            for (int iel = 0; iel < els_p4().size(); iel++) {
                pt = els_p4()[iel].pt();
                if (pt < 15.) continue;
                eta = els_p4()[iel].eta();
                if (fabs(eta) > 2.5) continue;
                phi = els_p4()[iel].phi();
                id = -11*els_charge()[iel];
                // ssid = isGoodLeptonNoIso(id,iel);
                // if (!ssid) continue;
                if (isMC) {
                    truth = is_truth_matched(id,iel);
                    mother_id = get_mother_id(id,iel);
                }
                miniiso = elMiniRelIsoCMS3_EA(iel, gconf.ea_version);
                ptrel = getPtRel(id, iel, true, ssWhichCorr);
                ptratio = pt/closestJet(els_p4()[iel], 0.4, 3.0, ssWhichCorr).pt();
                pucorr = rho * elEA03(iel, gconf.ea_version) * pow(getMiniDR(pt)/0.3,2);

                miniiso_ch = els_miniIso_ch()[iel] / pt;
                miniiso_nh = els_miniIso_nh()[iel] / pt;
                miniiso_em = els_miniIso_em()[iel] / pt;
                miniiso_raw = miniiso_ch + miniiso_nh + miniiso_em;
                miniiso_ncorr = miniiso_nh + miniiso_em - pucorr/pt;

                nmiss = els_exp_innerlayers()[iel];

                // pass_iso = (miniiso<0.09) && (ptratio>0.85 || ptrel>9.2);
                pass_id = isGoodLeptonNoIso(id,iel);
                pass_iso = isGoodLepton(id,iel);

                hovere = els_hOverE()[iel];

                out_tree->Fill();
            }


        }//event loop

        delete file;
    }//file loop

    out_file->cd();
    out_tree->Write();
    out_file->Close();

    return 0;

}

