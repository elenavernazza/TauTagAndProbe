#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <string.h>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <cstdlib>
#include <algorithm>

using namespace std;

bool CheckGoodMuon (bool in_MuIso, int in_MuID)
  {
    bool isPFMuon = (in_MuID & 0x2) >> 1;
    bool isGlobalMuon = (in_MuID & 0x4) >> 2;
    bool isTrackerMuon = (in_MuID & 0x8) >> 3;
    bool isLooseMuon = isPFMuon && (isGlobalMuon | isTrackerMuon);
    if (in_MuIso && isLooseMuon) return true;
    else return false;
  }

bool CheckGoodJet (int in_JetID, TString JetIDType)
  {
    bool looseJetID = in_JetID == 1;
    bool tightJetID = in_JetID == 2;
    bool tightLepVetoJetID = in_JetID == 3;
    if (JetIDType == "None") return true;
    if (JetIDType == "looseJetID") return looseJetID;
    if (JetIDType == "tightJetID") return tightJetID;
    if (JetIDType == "tightLepVetoJetID") return tightLepVetoJetID;
    else return false;
  }

// Check if a given event passes selections for the acceptance: offline + online + matching for all objects involved in VBF trigger
void ReadVBFEvent (TTree* inTree, UInt_t i_ev, vector<Float_t>* JetPt1_way1, vector<Float_t>* JetPt2_way1, vector<Float_t>* Mjj_way1, vector<Float_t>* JetPt1_way2, vector<Float_t>* JetPt2_way2, vector<Float_t>* Mjj_way2, TString JetIDType, bool JetSel30)
  {

    ULong64_t       in_EventNumber =  0;
    Int_t           in_RunNumber =  0;
    Int_t           in_lumi =  0;
    vector<float>   *in_PxJet = 0;
    vector<float>   *in_PyJet = 0;
    vector<float>   *in_PzJet = 0;
    vector<float>   *in_IDJet = 0;
    vector<float>   *in_l1tPtJet = 0;
    vector<float>   *in_l1tEtaJet = 0;
    vector<float>   *in_l1tPhiJet = 0;
    vector<float>   *in_MuPt = 0;
    vector<float>   *in_MuEta = 0;
    vector<float>   *in_MuPhi = 0;
    vector<bool>    *in_MuIso = 0;
    vector<int>     *in_MuID = 0;
    vector<float>   *in_l1tMuPt = 0;
    vector<float>   *in_l1tMuEta = 0;
    vector<float>   *in_l1tMuPhi = 0;
    vector<float>   *in_tauPt = 0;
    vector<float>   *in_tauEta = 0;
    vector<float>   *in_tauPhi = 0;
    vector<float>   *in_l1tauPt = 0;
    vector<float>   *in_l1tauEta = 0;
    vector<float>   *in_l1tauPhi = 0;

    inTree->SetBranchAddress("EventNumber", &in_EventNumber);
    inTree->SetBranchAddress("RunNumber", &in_RunNumber);
    inTree->SetBranchAddress("lumi", &in_lumi);
    inTree->SetBranchAddress("jets_px", &in_PxJet);
    inTree->SetBranchAddress("jets_py", &in_PyJet);
    inTree->SetBranchAddress("jets_pz", &in_PzJet);
    inTree->SetBranchAddress("jets_ID", &in_IDJet);
    inTree->SetBranchAddress("l1tPtJet", &in_l1tPtJet);
    inTree->SetBranchAddress("l1tEtaJet", &in_l1tEtaJet);
    inTree->SetBranchAddress("l1tPhiJet", &in_l1tPhiJet);
    inTree->SetBranchAddress("muons_pt", &in_MuPt);
    inTree->SetBranchAddress("muons_eta", &in_MuEta);
    inTree->SetBranchAddress("muons_phi", &in_MuPhi);
    inTree->SetBranchAddress("muons_PFIsoTight", &in_MuIso);
    inTree->SetBranchAddress("muons_type", &in_MuID);
    inTree->SetBranchAddress("l1t_muons_pt", &in_l1tMuPt);
    inTree->SetBranchAddress("l1t_muons_eta", &in_l1tMuEta);
    inTree->SetBranchAddress("l1t_muons_phi", &in_l1tMuPhi);
    inTree->SetBranchAddress("tauPt", &in_tauPt);
    inTree->SetBranchAddress("tauEta", &in_tauEta);
    inTree->SetBranchAddress("tauPhi", &in_tauPhi);
    inTree->SetBranchAddress("l1tPt", &in_l1tauPt);
    inTree->SetBranchAddress("l1tEta", &in_l1tauEta);
    inTree->SetBranchAddress("l1tPhi", &in_l1tauPhi);

    inTree->GetEntry(i_ev);

    // remove taus that are matched with jets within deltaR < 0.5
    vector <bool> jet_is_matched_to_tau;
    jet_is_matched_to_tau.clear();

    for (UInt_t ijet = 0; ijet < in_PxJet->size(); ++ijet)
      {
        TLorentzVector local_jet;
        Float_t jet_energy = sqrt(pow(in_PxJet->at(ijet),2)+pow(in_PyJet->at(ijet),2)+pow(in_PzJet->at(ijet),2));
        local_jet.SetPxPyPzE(in_PxJet->at(ijet), in_PyJet->at(ijet), in_PzJet->at(ijet), jet_energy);

        bool isMatched = false;
        for (UInt_t itau = 0; itau < in_tauPt->size(); ++itau)
          {
            TLorentzVector local_tau;
            local_tau.SetPtEtaPhiM(in_tauPt->at(itau), in_tauEta->at(itau), in_tauPhi->at(itau), 0.);
            if (local_tau.Pt()<10.) continue;
            if (local_tau.DeltaR(local_jet) < 0.5)
              {
                isMatched = true;
                break;
              }
          }
        jet_is_matched_to_tau.push_back(isMatched);
      }

    // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
    bool offline_objects_VBF = in_PxJet->size() > 1 && in_MuPt->size() > 0;
    bool offline_objects_MuTau = in_tauPt->size() > 0 ;
    bool L1_objects_VBF = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 ;

    if (offline_objects_VBF && L1_objects_VBF && offline_objects_MuTau)
      {

        bool check_VBF = false; // this has to be inside the loop
        bool check_VBF_online = false;
        bool check_VBF_offline = false;
        bool check_mu_matching = false;
        bool check_jet1_matching_way1 = false;
        bool check_jet2_matching_way1 = false;
        bool check_jet1_matching_way2 = false;
        bool check_jet2_matching_way2 = false;
        bool check_jets_VBF = false;

        TLorentzVector myGoodOfflineMuon;
        TLorentzVector myGoodOfflineJet1_way1;
        TLorentzVector myGoodOfflineJet2_way1;
        TLorentzVector myGoodOfflineDiJet_way1;
        TLorentzVector myGoodOfflineJet1_way2;
        TLorentzVector myGoodOfflineJet2_way2;
        TLorentzVector myGoodOfflineDiJet_way2;

        TLorentzVector myGoodOnlineMuon;
        TLorentzVector myGoodOnlineJet1_way1;
        TLorentzVector myGoodOnlineJet2_way1;
        TLorentzVector myGoodOnlineDiJet_way1;
        TLorentzVector myGoodOnlineJet1_way2;
        TLorentzVector myGoodOnlineJet2_way2;
        TLorentzVector myGoodOnlineDiJet_way2;

        TLorentzVector myGoodOfflineTau;
        myGoodOfflineTau.SetPtEtaPhiM(in_tauPt->at(0), in_tauEta->at(0), in_tauPhi->at(0), 1.776);

        for (UInt_t i_mu = 0 ; i_mu < in_MuPt->size() ; ++i_mu)
          {
            if (check_mu_matching) break;
            if (CheckGoodMuon(in_MuIso->at(i_mu), in_MuID->at(i_mu))) // Isolation and identification of muons
              {
                TLorentzVector myOfflineMuon;
                myOfflineMuon.SetPtEtaPhiM(in_MuPt->at(i_mu), in_MuEta->at(i_mu), in_MuPhi->at(i_mu), 0.105);

                for (UInt_t i_L1_mu = 0 ; i_L1_mu < in_l1tMuPt->size() ; ++i_L1_mu)
                  {
                    if (check_mu_matching) break;
                    TLorentzVector myOnlineMuon;
                    myOnlineMuon.SetPtEtaPhiM(in_l1tMuPt->at(i_L1_mu), in_l1tMuEta->at(i_L1_mu), in_l1tMuPhi->at(i_L1_mu), 0.105);
                    if (myOfflineMuon.DeltaR(myOnlineMuon) < 0.3)
                      {
                        myGoodOfflineMuon = myOfflineMuon;
                        myGoodOnlineMuon = myOnlineMuon;
                        check_mu_matching = true;
                        break;
                      }
                  }
              }
          }

        Float_t highest_mjj_offline_way1 = -99;
        int myGoodOfflineJet1Index_way1 = -1;
        int myGoodOfflineJet2Index_way1 = -1;
        int myGoodOnlineJet1Index_way1 = -1;
        int myGoodOnlineJet2Index_way1 = -1;
        Float_t highest_mjj_offline_way2 = -99;
        int myGoodOfflineJet1Index_way2 = -1;
        int myGoodOfflineJet2Index_way2 = -1;
        int myGoodOnlineJet1Index_way2 = -1;
        int myGoodOnlineJet2Index_way2 = -1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////// Way 1 /////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        // find the two jets giving highest mjj among jets with pt > 30
        for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size() ; ++i_jet1)
          {
            if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
            if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;
            Float_t jet1_energy_way1 = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
            TLorentzVector myOfflineJet1_way1;
            myOfflineJet1_way1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy_way1);
            if (JetSel30) {if (myOfflineJet1_way1.Pt() < 30) continue;}
            for (UInt_t i_jet2 = i_jet1 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
              {
                if (jet_is_matched_to_tau.at(i_jet2) == true) continue;
                if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
                Float_t jet2_energy_way1 = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                TLorentzVector myOfflineJet2_way1;
                myOfflineJet2_way1.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy_way1);
                if (JetSel30) {if (myOfflineJet2_way1.Pt() < 30) continue;}
                TLorentzVector myOfflineDiJet_way1;
                myOfflineDiJet_way1 = myOfflineJet1_way1 + myOfflineJet2_way1;
                if (myOfflineDiJet_way1.M() > highest_mjj_offline_way1)
                  {
                    myGoodOfflineJet1_way1 = myOfflineJet1_way1;
                    myGoodOfflineJet1Index_way1 = i_jet1;
                    myGoodOfflineJet2_way1 = myOfflineJet2_way1;
                    myGoodOfflineJet2Index_way1 = i_jet2;
                    myGoodOfflineDiJet_way1 = myGoodOfflineJet1_way1 + myGoodOfflineJet2_way1;
                    highest_mjj_offline_way1 = myGoodOfflineDiJet_way1.M();
                  }
              }
          }

        for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
          {
            if (check_jet1_matching_way1) break;
            TLorentzVector myOnlineJet1_way1;
            myOnlineJet1_way1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);
            if (myGoodOfflineJet1_way1.DeltaR(myOnlineJet1_way1) < 0.5)
              {
                myGoodOnlineJet1_way1 = myOnlineJet1_way1;
                myGoodOnlineJet1Index_way1 = i_L1_jet1;
                check_jet1_matching_way1 = true;
                break;
              }
          }

        for (UInt_t i_L1_jet2 = 0 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
          {
            if (check_jet2_matching_way1) break;
            if (int(i_L1_jet2) != myGoodOnlineJet1Index_way1)
              {
                TLorentzVector myOnlineJet2_way1;
                myOnlineJet2_way1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);
                if (myGoodOfflineJet2_way1.DeltaR(myOnlineJet2_way1) < 0.5)
                  {
                    myGoodOnlineJet2_way1 = myOnlineJet2_way1;
                    myGoodOnlineJet2Index_way1 = i_L1_jet2;
                    check_jet2_matching_way1 = true;
                    break;
                  }
              }
          }

        myGoodOnlineDiJet_way1 = myGoodOnlineJet1_way1 + myGoodOnlineJet2_way1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////// Way 2 /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size()-1 ; ++i_jet1)
          {
            if (check_jet1_matching_way2 && check_jet2_matching_way2) break;
            if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
            if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;

            // compute jet 1 energy from px, py, pz information assuming jet mass to 0
            Float_t jet1_energy_way2 = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
            // definition of Lorentz vector for offline jet 1 (px, py, pz, E [GeV])
            TLorentzVector myOfflineJet1_way2;
            myOfflineJet1_way2.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy_way2);
            if (JetSel30) {if (myOfflineJet1_way2.Pt() < 30) continue;}

            for (UInt_t i_jet2 = i_jet1+1 ; i_jet2 < in_PxJet->size() ; ++i_jet2)
              {
                if (check_jet1_matching_way2 && check_jet2_matching_way2) break;
                if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
                // compute jet 2 energy from px, py, pz information assuming jet mass to 0
                Float_t jet2_energy_way2 = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                // definition of Lorentz vector for offline jet 2 (px, py, pz, E [GeV])
                TLorentzVector myOfflineJet2_way2;
                myOfflineJet2_way2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy_way2);
                if (JetSel30) {if (myOfflineJet2_way2.Pt() < 30) continue;}

                for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
                  {

                    if (check_jet1_matching_way2 && check_jet2_matching_way2) break;
                    // definition of Lorentz vector for online jet 1 (pt, eta, phi, mass [GeV])
                    TLorentzVector myOnlineJet1_way2;
                    myOnlineJet1_way2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);

                    for (UInt_t i_L1_jet2 = i_L1_jet1+1 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
                      {

                        if (check_jet1_matching_way2 && check_jet2_matching_way2) break;
                        // definition of Lorentz vector for online jet 2 (pt, eta, phi, mass [GeV])
                        TLorentzVector myOnlineJet2_way2;
                        myOnlineJet2_way2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);

                        // check matching based on deltaR between online and offline object
                        if (myOfflineJet1_way2.DeltaR(myOnlineJet1_way2) < 0.5 && myOfflineJet2_way2.DeltaR(myOnlineJet2_way2) < 0.5)
                          {
                            myGoodOfflineJet1_way2 = myOfflineJet1_way2;
                            myGoodOfflineJet1Index_way2 = i_jet1;
                            myGoodOfflineJet2_way2 = myOfflineJet2_way2;
                            myGoodOfflineJet2Index_way2 = i_jet2;
                            myGoodOfflineDiJet_way2 = myGoodOfflineJet1_way2 + myGoodOfflineJet2_way2;

                            myGoodOnlineJet1_way2 = myOnlineJet1_way2;
                            myGoodOnlineJet1Index_way2 = i_L1_jet1;
                            myGoodOnlineJet2_way2 = myOnlineJet2_way2;
                            myGoodOnlineJet2Index_way2 = i_L1_jet2;
                            myGoodOnlineDiJet_way2 = myGoodOnlineJet1_way2 + myGoodOnlineJet2_way2;

                            check_jet1_matching_way2 = true;
                            check_jet2_matching_way2 = true;
                            break;
                          }
                      }
                  }
              }
          }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (check_jet1_matching_way1 && check_jet2_matching_way1)
          {
            JetPt1_way1->push_back(myGoodOfflineJet1_way1.Pt());
            JetPt2_way1->push_back(myGoodOfflineJet2_way1.Pt());
            Mjj_way1->push_back(myGoodOfflineDiJet_way1.M());
          }
        if (check_jet1_matching_way2 && check_jet2_matching_way2)
          {
            JetPt1_way2->push_back(myGoodOfflineJet1_way2.Pt());
            JetPt2_way2->push_back(myGoodOfflineJet2_way2.Pt());
            Mjj_way2->push_back(myGoodOfflineDiJet_way2.M());
          }
      }
  }

// Plot jetPt and mjj distribution for VBF events
void PlotDistributions(TString EventSample)
  {

    // TString EventSample = "MC_MiniAOD_VBFHHTo2B2Tau_12_01_23";
    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/"+EventSample;
    UInt_t events_Total = 0; 
    vector<Float_t> JetPt1_way1;
    vector<Float_t> JetPt2_way1;
    vector<Float_t> Mjj_way1;
    vector<Float_t> JetPt1_way2;
    vector<Float_t> JetPt2_way2;
    vector<Float_t> Mjj_way2;

    TString JetIDType = "tightLepVetoJetID";
    TString JetSel30 = "JetSel30";
    bool checkJetSel30 = false;
    if (JetSel30 == "JetSel30") {checkJetSel30 = true;}
    else if (JetSel30 == "NoJetSel30") {checkJetSel30 = false;}
    TString Output = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotDistributions/"+EventSample+"_2Ways_"+JetSel30+"_"+JetIDType;
    system("mkdir -p "+Output);

    // loop over all events in the ntuples to check which one passes the VBF or MuTau selections
    for (int nt = 0 ; nt < 5 ; ++nt)
      {
        TString FileName_ntuples = Path_ntuples + "/Ntuple_" + to_string(nt) + ".root";
        TFile f_in(FileName_ntuples.Data(), "READ");
        TDirectoryFile* df = (TDirectoryFile*)f_in.Get("Ntuplizer_noTagAndProbe_multipleTaus");
        TTree* inTree = (TTree*)df->Get("TagAndProbe");

        cout << "\nNumber of events in the Ntuple " << nt << " = " << inTree->GetEntries() << endl;
        events_Total += inTree->GetEntries();

        // each event is analysed for all the possible good combinations of selections
        for (UInt_t i_ev = 0 ; i_ev < inTree->GetEntries() ; ++i_ev)
          {
            ReadVBFEvent(inTree, i_ev, &JetPt1_way1, &JetPt2_way1, &Mjj_way1, &JetPt1_way2, &JetPt2_way2, &Mjj_way2, JetIDType, checkJetSel30);
          }
      }

    auto min_JetPt1_way1 = *min_element(JetPt1_way1.begin(), JetPt1_way1.end());
    auto max_JetPt1_way1 = *max_element(JetPt1_way1.begin(), JetPt1_way1.end());  
    TH1F* Histo_JetPt1_way1 = new TH1F("Histo_JetPt1_way1", "Histo_JetPt1_way1", 50, min_JetPt1_way1, max_JetPt1_way1);
    for (UInt_t i = 0; i < JetPt1_way1.size(); ++i)
      {
        Histo_JetPt1_way1->Fill(JetPt1_way1.at(i));
      }
 
    TH1F* Histo_JetPt2_way1 = new TH1F("Histo_JetPt2_way1", "Histo_JetPt2_way1", 50, min_JetPt1_way1, max_JetPt1_way1);
    for (UInt_t i = 0; i < JetPt2_way1.size(); ++i)
      {
        Histo_JetPt2_way1->Fill(JetPt2_way1.at(i));
      }

    auto min_JetPt1_way2 = *min_element(JetPt1_way2.begin(), JetPt1_way2.end());
    auto max_JetPt1_way2 = *max_element(JetPt1_way2.begin(), JetPt1_way2.end());  
    TH1F* Histo_JetPt1_way2 = new TH1F("Histo_JetPt1_way2", "Histo_JetPt1_way2", 50, min_JetPt1_way2, max_JetPt1_way2);
    for (UInt_t i = 0; i < JetPt1_way2.size(); ++i)
      {
        Histo_JetPt1_way2->Fill(JetPt1_way2.at(i));
      }

    TH1F* Histo_JetPt2_way2 = new TH1F("Histo_JetPt2_way2", "Histo_JetPt2_way2", 50, min_JetPt1_way2, max_JetPt1_way2);
    for (UInt_t i = 0; i < JetPt2_way2.size(); ++i)
      {
        Histo_JetPt2_way2->Fill(JetPt2_way2.at(i));
      }

    // gStyle->SetPalette(1,0);
    // auto cols = TColor::GetPalette();
    // cout << cols.GetSize() << endl;
    TCanvas c1 ("c1", "c1", 900., 650.);
    c1.SetGrid(10,10);
    Histo_JetPt1_way1->SetTitle("");
    Histo_JetPt1_way1->SetStats(0);

    Histo_JetPt1_way1->Draw();
    Histo_JetPt1_way1->SetLineWidth(3);
    Histo_JetPt1_way1->SetLineColor(kMagenta);

    Histo_JetPt1_way2->Draw("same");
    Histo_JetPt1_way2->SetLineWidth(3);
    Histo_JetPt1_way2->SetLineColor(kGreen+1);

    Histo_JetPt1_way1->GetXaxis()->SetTitle("Offline Jet p_{T} [GeV]");
    Histo_JetPt1_way1->GetXaxis()->SetRangeUser(0,500);
    // Histo_JetPt1_way1->GetXaxis()->SetTitleOffset(1.3);
    Histo_JetPt1_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_JetPt1_way1->GetYaxis()->SetTitle("entries");
    Histo_JetPt1_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_JetPt1_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_JetPt1_way1->GetYaxis()->SetRangeUser(0,1.05*Histo_JetPt1_way2->GetMaximum());

    TLegend* Legend =  new TLegend(0.65,0.72,0.88,0.88);
    Legend->SetBorderSize(0);
    Legend->AddEntry(Histo_JetPt1_way1 , "Way1: jet1", "LPE");
    Legend->AddEntry(Histo_JetPt1_way2 , "Way2: jet1", "LPE");
    Legend->Draw();

    TLatex Tex1;
    Tex1.SetTextSize(0.03);
    Tex1.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex1.Draw("same");

    TLatex Tex2;
    Tex2.SetTextSize(0.035);
    Tex2.SetTextAlign(31);
    Tex2.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex2.Draw("same");

    c1.SaveAs(Output+"/Distribution_JetPt1.png");
    c1.SaveAs(Output+"/Distribution_JetPt1.pdf");

    Histo_JetPt2_way1->Draw("same");
    Histo_JetPt2_way1->SetLineWidth(3);
    Histo_JetPt2_way1->SetLineColor(kMagenta+2);

    Histo_JetPt2_way2->Draw("same");
    Histo_JetPt2_way2->SetLineWidth(3);
    Histo_JetPt2_way2->SetLineColor(kGreen+3);
    Histo_JetPt1_way1->GetYaxis()->SetRangeUser(0,1.05*Histo_JetPt2_way1->GetMaximum());

    Legend->AddEntry(Histo_JetPt2_way1 , "Way1: jet2", "LPE");
    Legend->AddEntry(Histo_JetPt2_way2 , "Way2: jet2", "LPE");
    Legend->Draw();

    c1.SaveAs(Output+"/Distribution_JetPt.png");
    c1.SaveAs(Output+"/Distribution_JetPt.pdf");
    c1.Close();

    TCanvas c2 ("c2", "c2", 900., 650.);
    c2.SetGrid(10,10);
    Histo_JetPt2_way1->SetTitle("");
    Histo_JetPt2_way1->SetStats(0);

    Histo_JetPt2_way1->Draw();
    Histo_JetPt2_way1->SetLineWidth(3);
    Histo_JetPt2_way1->SetLineColor(kMagenta);

    Histo_JetPt2_way2->Draw("same");
    Histo_JetPt2_way2->SetLineWidth(3);
    Histo_JetPt2_way2->SetLineColor(kGreen+1);

    Histo_JetPt2_way1->GetXaxis()->SetTitle("Offline Jet p_{T} [GeV]");
    Histo_JetPt2_way1->GetXaxis()->SetRangeUser(0,500);
    // Histo_JetPt2_way1->GetXaxis()->SetTitleOffset(1.3);
    Histo_JetPt2_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_JetPt2_way1->GetYaxis()->SetTitle("entries");
    Histo_JetPt2_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_JetPt2_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_JetPt2_way1->GetYaxis()->SetRangeUser(0,1.05*Histo_JetPt2_way1->GetMaximum());

    TLegend* Legend1 =  new TLegend(0.65,0.72,0.88,0.88);
    Legend1->SetBorderSize(0);
    Legend1->AddEntry(Histo_JetPt2_way1 , "Way1: jet2", "LPE");
    Legend1->AddEntry(Histo_JetPt2_way2 , "Way2: jet2", "LPE");
    Legend1->Draw();

    TLatex Tex11;
    Tex11.SetTextSize(0.03);
    Tex11.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex11.Draw("same");

    TLatex Tex22;
    Tex22.SetTextSize(0.035);
    Tex22.SetTextAlign(31);
    Tex22.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex22.Draw("same");

    c2.SaveAs(Output+"/Distribution_JetPt2.png");
    c2.SaveAs(Output+"/Distribution_JetPt2.pdf");
    c2.Close();

    // plot distributions Mjj
    auto min_Mjj_way1 = *min_element(Mjj_way1.begin(), Mjj_way1.end());
    auto max_Mjj_way1 = *max_element(Mjj_way1.begin(), Mjj_way1.end());  
    TH1F* Histo_Mjj_way1 = new TH1F("Histo_Mjj_way1", "Histo_Mjj_way1", 50, min_Mjj_way1, max_Mjj_way1);
    for (UInt_t i = 0; i < Mjj_way1.size(); ++i)
      {
        Histo_Mjj_way1->Fill(Mjj_way1.at(i));
      }

    // plot distributions Mjj
    auto min_Mjj_way2 = *min_element(Mjj_way2.begin(), Mjj_way2.end());
    auto max_Mjj_way2 = *max_element(Mjj_way2.begin(), Mjj_way2.end());  
    TH1F* Histo_Mjj_way2 = new TH1F("Histo_Mjj_way2", "Histo_Mjj_way2", 50, min_Mjj_way2, max_Mjj_way2);
    for (UInt_t i = 0; i < Mjj_way2.size(); ++i)
      {
        Histo_Mjj_way2->Fill(Mjj_way2.at(i));
      }

    TCanvas c3 ("c3", "c3", 900., 650.);
    c3.SetGrid(10,10);
    Histo_Mjj_way1->SetTitle("");
    Histo_Mjj_way1->SetStats(0);

    Histo_Mjj_way1->Draw();
    Histo_Mjj_way1->SetLineWidth(3);
    Histo_Mjj_way1->SetLineColor(kMagenta);

    Histo_Mjj_way2->Draw("same");
    Histo_Mjj_way2->SetLineWidth(3);
    Histo_Mjj_way2->SetLineColor(kGreen+1);

    Histo_Mjj_way1->GetXaxis()->SetTitle("M_{jj} [GeV]");
    // Histo_Mjj_way1->GetXaxis()->SetTitleOffset(1.3);
    Histo_Mjj_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_Mjj_way1->GetYaxis()->SetTitle("entries");
    Histo_Mjj_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_Mjj_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_Mjj_way1->GetYaxis()->SetRangeUser(0,1.05*Histo_Mjj_way2->GetMaximum());
    Histo_Mjj_way1->GetXaxis()->SetRangeUser(0,5000);

    TLegend* Legend2 =  new TLegend(0.65,0.8,0.88,0.88);
    Legend2->SetBorderSize(0);
    Legend2->AddEntry(Histo_Mjj_way1 , "Way1: mjj", "LPE");
    Legend2->AddEntry(Histo_Mjj_way2 , "Way2: mjj", "LPE");
    Legend2->Draw();

    TLatex Tex3;
    Tex3.SetTextSize(0.03);
    Tex3.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex3.Draw("same");

    TLatex Tex4;
    Tex4.SetTextSize(0.035);
    Tex4.SetTextAlign(31);
    Tex4.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex4.Draw("same");

    c3.SaveAs(Output+"/Distribution_Mjj.png");
    c3.SaveAs(Output+"/Distribution_Mjj.pdf");
    c3.Close();

  }