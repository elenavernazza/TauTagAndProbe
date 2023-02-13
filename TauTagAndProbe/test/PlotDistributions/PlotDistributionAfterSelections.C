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
void ReadVBFEvent (TTree* inTree, UInt_t i_ev, vector<Float_t>* JetPt1, vector<Float_t>* JetPt2, vector<Float_t>* Mjj, vector<Float_t>* MuPt, vector<Float_t>* TauPt, TString JetIDType, bool JetSel30, bool checkmu, bool checkjets, TString Method)
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
        bool check_jet1_matching = false;
        bool check_jet2_matching = false;

        TLorentzVector myGoodOfflineMuon;
        TLorentzVector myGoodOfflineJet1;
        TLorentzVector myGoodOfflineJet2;
        TLorentzVector myGoodOfflineDiJet;

        TLorentzVector myGoodOnlineMuon;
        TLorentzVector myGoodOnlineJet1;
        TLorentzVector myGoodOnlineJet2;
        TLorentzVector myGoodOnlineDiJet;

        TLorentzVector myGoodOfflineTau;
        myGoodOfflineTau.SetPtEtaPhiM(in_tauPt->at(0), in_tauEta->at(0), in_tauPhi->at(0), 1.776);

        for (UInt_t i_mu = 0 ; i_mu < in_MuPt->size() ; ++i_mu)
          {
            if (check_mu_matching) break;
            if (checkmu) {if (!CheckGoodMuon(in_MuIso->at(i_mu), in_MuID->at(i_mu))) continue;}
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

        Float_t highest_mjj_offline = -99;
        int myGoodOfflineJet1Index = -1;
        int myGoodOfflineJet2Index = -1;
        int myGoodOnlineJet1Index = -1;
        int myGoodOnlineJet2Index = -1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////// Way 1 /////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        if (Method == "Way1")
          {
            // find the two jets giving highest mjj among jets with pt > 30
            for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size() ; ++i_jet1)
              {
                if (checkjets)
                  {
                    if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
                    if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;
                  }
                Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                TLorentzVector myOfflineJet1;
                myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
                if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}
                for (UInt_t i_jet2 = i_jet1 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
                  {
                    if (checkjets)
                      {
                        if (jet_is_matched_to_tau.at(i_jet2) == true) continue;
                        if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
                      }
                    Float_t jet2_energy = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                    TLorentzVector myOfflineJet2;
                    myOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy);
                    if (JetSel30) {if (myOfflineJet2.Pt() < 30) continue;}
                    TLorentzVector myOfflineDiJet;
                    myOfflineDiJet = myOfflineJet1 + myOfflineJet2;
                    if (myOfflineDiJet.M() > highest_mjj_offline)
                      {
                        myGoodOfflineJet1 = myOfflineJet1;
                        myGoodOfflineJet1Index = i_jet1;
                        myGoodOfflineJet2 = myOfflineJet2;
                        myGoodOfflineJet2Index = i_jet2;
                        myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;
                        highest_mjj_offline = myGoodOfflineDiJet.M();
                      }
                  }
              }

            for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
              {
                if (check_jet1_matching) break;
                TLorentzVector myOnlineJet1;
                myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);
                if (myGoodOfflineJet1.DeltaR(myOnlineJet1) < 0.5)
                  {
                    myGoodOnlineJet1 = myOnlineJet1;
                    myGoodOnlineJet1Index = i_L1_jet1;
                    check_jet1_matching = true;
                    break;
                  }
              }

            for (UInt_t i_L1_jet2 = 0 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
              {
                if (check_jet2_matching) break;
                if (int(i_L1_jet2) != myGoodOnlineJet1Index)
                  {
                    TLorentzVector myOnlineJet2;
                    myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);
                    if (myGoodOfflineJet2.DeltaR(myOnlineJet2) < 0.5)
                      {
                        myGoodOnlineJet2 = myOnlineJet2;
                        myGoodOnlineJet2Index = i_L1_jet2;
                        check_jet2_matching = true;
                        break;
                      }
                  }
              }

            myGoodOnlineDiJet = myGoodOnlineJet1 + myGoodOnlineJet2;

          }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////// Way 2 /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (Method == "Way2")
          {

            for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size()-1 ; ++i_jet1)
              {
                if (check_jet1_matching && check_jet2_matching) break;
                if (checkjets)
                  {
                    if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
                    if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;
                  }
                // compute jet 1 energy from px, py, pz information assuming jet mass to 0
                Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                // definition of Lorentz vector for offline jet 1 (px, py, pz, E [GeV])
                TLorentzVector myOfflineJet1;
                myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
                if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}

                for (UInt_t i_jet2 = i_jet1+1 ; i_jet2 < in_PxJet->size() ; ++i_jet2)
                  {
                    if (checkjets)
                      {
                        if (check_jet1_matching && check_jet2_matching) break;
                        if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
                      }
                    // compute jet 2 energy from px, py, pz information assuming jet mass to 0
                    Float_t jet2_energy = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                    // definition of Lorentz vector for offline jet 2 (px, py, pz, E [GeV])
                    TLorentzVector myOfflineJet2;
                    myOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy);
                    if (JetSel30) {if (myOfflineJet2.Pt() < 30) continue;}

                    for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
                      {

                        if (check_jet1_matching && check_jet2_matching) break;
                        // definition of Lorentz vector for online jet 1 (pt, eta, phi, mass [GeV])
                        TLorentzVector myOnlineJet1;
                        myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);

                        for (UInt_t i_L1_jet2 = i_L1_jet1+1 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
                          {

                            if (check_jet1_matching && check_jet2_matching) break;
                            // definition of Lorentz vector for online jet 2 (pt, eta, phi, mass [GeV])
                            TLorentzVector myOnlineJet2;
                            myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);

                            // check matching based on deltaR between online and offline object
                            if (myOfflineJet1.DeltaR(myOnlineJet1) < 0.5 && myOfflineJet2.DeltaR(myOnlineJet2) < 0.5)
                              {
                                myGoodOfflineJet1 = myOfflineJet1;
                                myGoodOfflineJet1Index = i_jet1;
                                myGoodOfflineJet2 = myOfflineJet2;
                                myGoodOfflineJet2Index = i_jet2;
                                myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;

                                myGoodOnlineJet1 = myOnlineJet1;
                                myGoodOnlineJet1Index = i_L1_jet1;
                                myGoodOnlineJet2 = myOnlineJet2;
                                myGoodOnlineJet2Index = i_L1_jet2;
                                myGoodOnlineDiJet = myGoodOnlineJet1 + myGoodOnlineJet2;

                                check_jet1_matching = true;
                                check_jet2_matching = true;
                                break;
                              }
                          }
                      }
                  }
              }
          }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (check_jet1_matching && check_jet2_matching)
          {
            JetPt1->push_back(myGoodOfflineJet1.Pt());
            JetPt2->push_back(myGoodOfflineJet2.Pt());
            Mjj->push_back(myGoodOfflineDiJet.M());
          }
        if (check_mu_matching)
          {
            MuPt->push_back(myGoodOfflineMuon.Pt());
          }
        if (myGoodOfflineTau.Pt() > 0)
          {
            TauPt->push_back(myGoodOfflineTau.Pt());
          }
      }
  }

// Plot jetPt and mjj distribution for VBF events
void PlotDistributionsAfterSelections(TString EventSample)
  {

    // TString EventSample = "MC_MiniAOD_VBFHHTo2B2Tau_12_01_23";
    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/"+EventSample;
    UInt_t events_Total = 0; 
    vector<Float_t> JetPt1_way1;
    vector<Float_t> JetPt2_way1;
    vector<Float_t> Mjj_way1;
    vector<Float_t> MuPt_way1;
    vector<Float_t> TauPt_way1;
    vector<Float_t> JetPt1_way2;
    vector<Float_t> JetPt2_way2;
    vector<Float_t> Mjj_way2;
    vector<Float_t> MuPt_way2;
    vector<Float_t> TauPt_way2;

    TString Method = "Way1";
    TString JetIDType = "tightLepVetoJetID";
    TString JetSel30 = "JetSel30";
    bool checkJetSel30 = false;
    if (JetSel30 == "JetSel30") {checkJetSel30 = true;}
    else if (JetSel30 == "NoJetSel30") {checkJetSel30 = false;}
    TString Output = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotDistributions/"+EventSample+"_DistributionAfterSelection_"+JetSel30+"_"+JetIDType;
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
            ReadVBFEvent(inTree, i_ev, &JetPt1_way1, &JetPt2_way1, &Mjj_way1, &MuPt_way1, &TauPt_way1, JetIDType, checkJetSel30, true, true, Method);
          }
      }

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
            ReadVBFEvent(inTree, i_ev, &JetPt1_way2, &JetPt2_way2, &Mjj_way2, &MuPt_way2, &TauPt_way2, JetIDType, checkJetSel30, false, false, Method);
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
  
    TH1F* Histo_JetPt1_way2 = new TH1F("Histo_JetPt1_way2", "Histo_JetPt1_way2", 50, min_JetPt1_way1, max_JetPt1_way1);
    for (UInt_t i = 0; i < JetPt1_way2.size(); ++i)
      {
        Histo_JetPt1_way2->Fill(JetPt1_way2.at(i));
      }

    TH1F* Histo_JetPt2_way2 = new TH1F("Histo_JetPt2_way2", "Histo_JetPt2_way2", 50, min_JetPt1_way1, max_JetPt1_way1);
    for (UInt_t i = 0; i < JetPt2_way2.size(); ++i)
      {
        Histo_JetPt2_way2->Fill(JetPt2_way2.at(i));
      }

    auto min_Mjj_way1 = *min_element(Mjj_way1.begin(), Mjj_way1.end());
    auto max_Mjj_way1 = *max_element(Mjj_way1.begin(), Mjj_way1.end());  
    TH1F* Histo_Mjj_way1 = new TH1F("Histo_Mjj_way1", "Histo_Mjj_way1", 50, min_Mjj_way1, max_Mjj_way1);
    for (UInt_t i = 0; i < Mjj_way1.size(); ++i)
      {
        Histo_Mjj_way1->Fill(Mjj_way1.at(i));
      }

    TH1F* Histo_Mjj_way2 = new TH1F("Histo_Mjj_way2", "Histo_Mjj_way2", 50, min_Mjj_way1, max_Mjj_way1);
    for (UInt_t i = 0; i < Mjj_way2.size(); ++i)
      {
        Histo_Mjj_way2->Fill(Mjj_way2.at(i));
      }

    auto min_MuPt_way1 = *min_element(MuPt_way1.begin(), MuPt_way1.end());
    auto max_MuPt_way1 = *max_element(MuPt_way1.begin(), MuPt_way1.end());  
    TH1F* Histo_MuPt_way1 = new TH1F("Histo_MuPt_way1", "Histo_MuPt_way1", 50, min_MuPt_way1, max_MuPt_way1);
    for (UInt_t i = 0; i < MuPt_way1.size(); ++i)
      {
        Histo_MuPt_way1->Fill(MuPt_way1.at(i));
      }

    TH1F* Histo_MuPt_way2 = new TH1F("Histo_MuPt_way2", "Histo_MuPt_way2", 50, min_MuPt_way1, max_MuPt_way1);
    for (UInt_t i = 0; i < MuPt_way2.size(); ++i)
      {
        Histo_MuPt_way2->Fill(MuPt_way2.at(i));
      }

    auto min_TauPt_way1 = *min_element(TauPt_way1.begin(), TauPt_way1.end());
    auto max_TauPt_way1 = *max_element(TauPt_way1.begin(), TauPt_way1.end());  
    TH1F* Histo_TauPt_way1 = new TH1F("Histo_TauPt_way1", "Histo_TauPt_way1", 50, min_TauPt_way1, max_TauPt_way1);
    for (UInt_t i = 0; i < TauPt_way1.size(); ++i)
      {
        Histo_TauPt_way1->Fill(TauPt_way1.at(i));
      }

    TH1F* Histo_TauPt_way2 = new TH1F("Histo_TauPt_way2", "Histo_TauPt_way2", 50, min_TauPt_way1, max_TauPt_way1);
    for (UInt_t i = 0; i < TauPt_way2.size(); ++i)
      {
        Histo_TauPt_way2->Fill(TauPt_way2.at(i));
      }

    ////////////////////////////////////////////////////////////////////////////////

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

    Histo_JetPt1_way1->GetXaxis()->SetTitle("Offline Jet p_{T}^{1} [GeV]");
    Histo_JetPt1_way1->GetXaxis()->SetRangeUser(0,500);
    Histo_JetPt1_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_JetPt1_way1->GetYaxis()->SetTitle("entries");
    Histo_JetPt1_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_JetPt1_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_JetPt1_way1->GetYaxis()->SetRangeUser(0,1.05*max(Histo_JetPt1_way1->GetMaximum(), Histo_JetPt1_way2->GetMaximum()));

    TLegend* Legend =  new TLegend(0.55,0.8,0.88,0.88);
    Legend->SetBorderSize(0);
    Legend->AddEntry(Histo_JetPt1_way1 , "No selections", "LPE");
    Legend->AddEntry(Histo_JetPt1_way2 , "TightJetID & TauVeto & LepVeto", "LPE");
    Legend->Draw();

    TLatex Tex1;
    Tex1.SetTextSize(0.03);
    Tex1.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex1.Draw("same");

    TLatex Tex11;
    Tex11.SetTextSize(0.035);
    Tex11.SetTextAlign(31);
    Tex11.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex11.Draw("same");

    c1.SaveAs(Output+"/Distribution_JetPt1.png");
    c1.SaveAs(Output+"/Distribution_JetPt1.pdf");
    c1.Close();

    ////////////////////////////////////////////////////////////////////////////////

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

    Histo_JetPt2_way1->GetXaxis()->SetTitle("Offline Jet p_{T}^{1} [GeV]");
    Histo_JetPt2_way1->GetXaxis()->SetRangeUser(0,500);
    Histo_JetPt2_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_JetPt2_way1->GetYaxis()->SetTitle("entries");
    Histo_JetPt2_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_JetPt2_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_JetPt2_way1->GetYaxis()->SetRangeUser(0,1.05*max(Histo_JetPt2_way1->GetMaximum(), Histo_JetPt2_way2->GetMaximum()));

    TLegend* Legend2 =  new TLegend(0.55,0.8,0.88,0.88);
    Legend2->SetBorderSize(0);
    Legend2->AddEntry(Histo_JetPt2_way1 , "No selections", "LPE");
    Legend2->AddEntry(Histo_JetPt2_way2 , "TightJetID & TauVeto & LepVeto", "LPE");
    Legend2->Draw();

    TLatex Tex2;
    Tex2.SetTextSize(0.03);
    Tex2.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex2.Draw("same");

    TLatex Tex22;
    Tex22.SetTextSize(0.035);
    Tex22.SetTextAlign(31);
    Tex22.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex22.Draw("same");

    c2.SaveAs(Output+"/Distribution_JetPt2.png");
    c2.SaveAs(Output+"/Distribution_JetPt2.pdf");
    c2.Close();

    ////////////////////////////////////////////////////////////////////////////////

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
    Histo_Mjj_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_Mjj_way1->GetYaxis()->SetTitle("entries");
    Histo_Mjj_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_Mjj_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_Mjj_way1->GetYaxis()->SetRangeUser(0,1.05*max(Histo_Mjj_way2->GetMaximum(), Histo_Mjj_way1->GetMaximum()));
    Histo_Mjj_way1->GetXaxis()->SetRangeUser(0,5000);

    TLegend* Legend3 =  new TLegend(0.55,0.8,0.88,0.88);
    Legend3->SetBorderSize(0);
    Legend3->AddEntry(Histo_Mjj_way1 , "No selections", "LPE");
    Legend3->AddEntry(Histo_Mjj_way2 , "TightJetID & TauVeto & LepVeto", "LPE");
    Legend3->Draw();

    TLatex Tex3;
    Tex3.SetTextSize(0.03);
    Tex3.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex3.Draw("same");

    TLatex Tex33;
    Tex33.SetTextSize(0.035);
    Tex33.SetTextAlign(31);
    Tex33.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex33.Draw("same");

    c3.SaveAs(Output+"/Distribution_Mjj.png");
    c3.SaveAs(Output+"/Distribution_Mjj.pdf");
    c3.Close();

    ////////////////////////////////////////////////////////////////////////////////

    TCanvas c4 ("c4", "c4", 900., 650.);
    c4.SetGrid(10,10);
    Histo_MuPt_way1->SetTitle("");
    Histo_MuPt_way1->SetStats(0);

    Histo_MuPt_way1->Draw();
    Histo_MuPt_way1->SetLineWidth(3);
    Histo_MuPt_way1->SetLineColor(kMagenta);

    Histo_MuPt_way2->Draw("same");
    Histo_MuPt_way2->SetLineWidth(3);
    Histo_MuPt_way2->SetLineColor(kGreen+1);

    Histo_MuPt_way1->GetXaxis()->SetTitle("Offline #mu p_{T} [GeV]");
    Histo_MuPt_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_MuPt_way1->GetYaxis()->SetTitle("entries");
    Histo_MuPt_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_MuPt_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_MuPt_way1->GetYaxis()->SetRangeUser(0,1.05*max(Histo_MuPt_way2->GetMaximum(), Histo_MuPt_way1->GetMaximum()));
    Histo_MuPt_way1->GetXaxis()->SetRangeUser(0,5000);

    TLegend* Legend4 =  new TLegend(0.55,0.8,0.88,0.88);
    Legend4->SetBorderSize(0);
    Legend4->AddEntry(Histo_MuPt_way1 , "No selections", "LPE");
    Legend4->AddEntry(Histo_MuPt_way2 , "Isolation & Identification", "LPE");
    Legend4->Draw();

    TLatex Tex4;
    Tex4.SetTextSize(0.03);
    Tex4.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex4.Draw("same");

    TLatex Tex44;
    Tex44.SetTextSize(0.035);
    Tex44.SetTextAlign(31);
    Tex44.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex44.Draw("same");

    c4.SaveAs(Output+"/Distribution_MuPt.png");
    c4.SaveAs(Output+"/Distribution_MuPt.pdf");
    c4.Close();

    ////////////////////////////////////////////////////////////////////////////////

    TCanvas c5 ("c5", "c5", 900., 650.);
    c5.SetGrid(10,10);
    Histo_TauPt_way1->SetTitle("");
    Histo_TauPt_way1->SetStats(0);

    Histo_TauPt_way1->Draw();
    Histo_TauPt_way1->SetLineWidth(3);
    Histo_TauPt_way1->SetLineColor(kMagenta);

    // Histo_TauPt_way2->Draw("same");
    // Histo_TauPt_way2->SetLineWidth(3);
    // Histo_TauPt_way2->SetLineColor(kGreen+1);

    Histo_TauPt_way1->GetXaxis()->SetTitle("Offline #tau p_{T} [GeV]");
    Histo_TauPt_way1->GetXaxis()->SetTitleSize(0.04);
    Histo_TauPt_way1->GetYaxis()->SetTitle("entries");
    Histo_TauPt_way1->GetYaxis()->SetTitleOffset(1.3);
    Histo_TauPt_way1->GetYaxis()->SetTitleSize(0.04);
    Histo_TauPt_way1->GetYaxis()->SetRangeUser(0,1.05*max(Histo_TauPt_way2->GetMaximum(), Histo_TauPt_way1->GetMaximum()));
    Histo_TauPt_way1->GetXaxis()->SetRangeUser(0,5000);

    TLegend* Legend5 =  new TLegend(0.55,0.84,0.88,0.88);
    Legend5->SetBorderSize(0);
    Legend5->AddEntry(Histo_TauPt_way1 , "No selections", "LPE");
    // Legend5->AddEntry(Histo_TauPt_way2 , "Isolation & Identification", "LPE");
    Legend5->Draw();

    TLatex Tex5;
    Tex5.SetTextSize(0.03);
    Tex5.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex5.Draw("same");

    TLatex Tex55;
    Tex55.SetTextSize(0.035);
    Tex55.SetTextAlign(31);
    Tex55.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex55.Draw("same");

    c5.SaveAs(Output+"/Distribution_TauPt.png");
    c5.SaveAs(Output+"/Distribution_TauPt.pdf");
    c5.Close();


  }