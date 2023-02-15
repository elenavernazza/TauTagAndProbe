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
#include <TPad.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <string.h>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <cstdlib>
#include <algorithm>
#include <TMultiGraph.h>
#include "../Utils/CheckL1Triggers.h"

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

bool CheckMuonQuality (int in_l1tMuQual)
  {
    // open quality: 4, 5, 6, 7, 8, 9, 10, 11, 12,13, 14, 15
    // single quality: 12, 13, 14, 15
    if (in_l1tMuQual > 3) { return true; }
    else { return false; }
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
void CheckVBF (TTree* inTree, UInt_t i_ev, vector<array<Float_t, 4>> set_of_on_cuts, vector<array<Float_t, 4>> set_of_off_cuts, TH1F* Total_Muon, TH1F* Match_Muon, TH1F* Offline_Muon, TH1F* L1_Offline_Muon, TH1F* Total_Jet1, TH1F* Match_Jet1, TH1F* Offline_Jet1, TH1F* L1_Offline_Jet1, TH1F* Total_Jet2, TH1F* Match_Jet2, TH1F* Offline_Jet2, TH1F* L1_Offline_Jet2, TH1F* Total_Mjj, TH1F* Match_Mjj, TH1F* Offline_Mjj, TH1F* L1_Offline_Mjj, bool JetSel30, TString JetIDType, TString Method)
  {

    bool checkmu = true;
    bool checkjets = true;

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
    vector<float>   *in_l1tMuQual = 0;
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
    inTree->SetBranchAddress("l1t_muons_qual", &in_l1tMuQual);
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
        for (UInt_t i_cut = 0; i_cut < set_of_on_cuts.size(); ++i_cut)
          {
            bool check_VBF = false;
            bool check_mu_matching = false;
            bool check_jet1_matching = false;
            bool check_jet2_matching = false;

            bool check_jet1_online = false;
            bool check_jet2_online = false;
            bool check_mjj_online = false;
            bool check_mu_online = false;
            
            bool check_jet1_offline = false;
            bool check_jet2_offline = false;
            bool check_mjj_offline = false;
            bool check_mu_offline = false;

            bool check_VBF_online = false;
            bool check_VBF_offline = false;

            bool check_VBF_offline_ptj1 = false;
            bool check_VBF_offline_ptj2 = false;
            bool check_VBF_offline_mjj = false;
            bool check_VBF_offline_muon = false;

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
                    if (CheckMuonQuality(int(in_l1tMuQual->at(i_L1_mu))) == false) continue;
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

            if (check_mu_matching && check_jet1_matching && check_jet2_matching)
              {
                check_jet1_online = myGoodOnlineJet1.Pt() > set_of_on_cuts.at(i_cut)[0];
                check_jet2_online = myGoodOnlineJet2.Pt() > set_of_on_cuts.at(i_cut)[1];
                check_mjj_online  = myGoodOnlineDiJet.M() > set_of_on_cuts.at(i_cut)[2];
                check_mu_online   = myGoodOnlineMuon.Pt() > set_of_on_cuts.at(i_cut)[3];
                
                check_jet1_offline = myGoodOfflineJet1.Pt() > set_of_off_cuts.at(i_cut)[0];
                check_jet2_offline = myGoodOfflineJet2.Pt() > set_of_off_cuts.at(i_cut)[1];
                check_mjj_offline  = myGoodOfflineDiJet.M() > set_of_off_cuts.at(i_cut)[2];
                check_mu_offline   = myGoodOfflineMuon.Pt() > set_of_off_cuts.at(i_cut)[3];

                check_VBF_online = check_jet1_online && check_jet2_online && check_mjj_online && check_mu_online;
                check_VBF_offline = check_jet1_offline && check_jet2_offline && check_mjj_offline && check_mu_offline;

                check_VBF_offline_ptj1 = check_jet2_offline && check_mjj_offline && check_mu_offline;
                check_VBF_offline_ptj2 = check_jet1_offline && check_mjj_offline && check_mu_offline;
                check_VBF_offline_mjj = check_jet1_offline && check_jet2_offline && check_mu_offline;
                check_VBF_offline_muon = check_jet1_offline && check_jet2_offline && check_mjj_offline;

                check_VBF = check_VBF_online && check_VBF_offline;
              }

            if (myGoodOfflineMuon.Pt() > 0.001)
              {
                Total_Muon->Fill(myGoodOfflineMuon.Pt());
                if (check_mu_matching) Match_Muon->Fill(myGoodOfflineMuon.Pt());
                if (check_VBF_offline_muon) Offline_Muon->Fill(myGoodOfflineMuon.Pt());
                if (check_VBF_online && check_VBF_offline_muon) L1_Offline_Muon->Fill(myGoodOfflineMuon.Pt());
              }

            if (myGoodOfflineJet1.Pt() > 0.001)
              {
                Total_Jet1->Fill(myGoodOfflineJet1.Pt());
                if (check_jet1_matching) Match_Jet1->Fill(myGoodOfflineJet1.Pt());
                if (check_VBF_offline_ptj1) Offline_Jet1->Fill(myGoodOfflineJet1.Pt());
                if (check_VBF_online && check_VBF_offline_ptj1) L1_Offline_Jet1->Fill(myGoodOfflineJet1.Pt());
              }

            if (myGoodOfflineJet2.Pt() > 0.001)
              {
                Total_Jet2->Fill(myGoodOfflineJet2.Pt());
                if (check_jet2_matching) Match_Jet2->Fill(myGoodOfflineJet2.Pt());
                if (check_VBF_offline_ptj2) Offline_Jet2->Fill(myGoodOfflineJet2.Pt());
                if (check_VBF_online && check_VBF_offline_ptj2) L1_Offline_Jet2->Fill(myGoodOfflineJet2.Pt());
              }

            if (myGoodOfflineDiJet.M() > 0.0001)
              {
                Total_Mjj->Fill(myGoodOfflineDiJet.M());
                if (check_jet1_matching && check_jet2_matching) Match_Mjj->Fill(myGoodOfflineDiJet.M());
                if (check_VBF_offline_mjj) Offline_Mjj->Fill(myGoodOfflineDiJet.M());
                if (check_VBF_online && check_VBF_offline_mjj) L1_Offline_Mjj->Fill(myGoodOfflineDiJet.M());
              }
        }
      }
  }

void PlotDistributions (TString output, TH1F* Total, TH1F* Match, TH1F* Offline, TH1F* L1_Offline, TString DistributionName, TString XaxisTitle)
  {
    
    TString rootname = output+"/"+DistributionName+".root";
    TFile* out_file = TFile::Open(rootname, "recreate");
    out_file->cd();

    TCanvas c ("c", "c", 900., 900.);
    c.cd();
    c.SetGrid(10,10);

    Total->SetTitle("");
    Total->GetYaxis()->SetRangeUser(0,1.05*Total->GetBinContent(Total->GetMaximumBin()));
    Match->GetYaxis()->SetRangeUser(0,1.05*Total->GetBinContent(Total->GetMaximumBin()));
    L1_Offline->GetYaxis()->SetRangeUser(0,1.05*Total->GetBinContent(Total->GetMaximumBin()));
    Offline->GetYaxis()->SetRangeUser(0,1.05*Total->GetBinContent(Total->GetMaximumBin()));
    Total->SetStats(0);
    Int_t steps = (Total->GetXaxis()->GetXmax() - Total->GetXaxis()->GetXmin())/Total->GetNbinsX();
    TString title = "Entries/"+to_string(steps)+" GeV";
    Total->GetYaxis()->SetTitle(title);
    Total->GetYaxis()->SetTitleOffset(1.55);
    Total->GetXaxis()->SetTitle(XaxisTitle);
    Total->GetXaxis()->SetTitleOffset(1.2);
    Total->Draw();
    Total->SetLineWidth(3);
    Total->SetLineColor(kAzure);
    Total->SetLineStyle(7);
    Match->Draw("same");
    Match->SetLineWidth(3);
    Match->SetLineStyle(1);
    Match->SetLineColor(kAzure+5);
    Match->SetFillColorAlpha(kAzure+5, 0.05);
    Offline->Draw("same");
    Offline->SetLineWidth(3);
    Offline->SetLineStyle(7);
    Offline->SetLineColor(kViolet);
    L1_Offline->Draw("same");
    L1_Offline->SetLineWidth(3);
    L1_Offline->SetLineStyle(1);
    L1_Offline->SetLineColor(kViolet+5);
    L1_Offline->SetFillColorAlpha(kViolet+5, 0.05);

    Match->Write();
    Total->Write();
    L1_Offline->Write();
    Offline->Write();

    TLegend Legend11 (0.45,0.74,0.89,0.89);
    Legend11.SetBorderSize(0);
    Legend11.AddEntry(Total, "Inclusive", "LPE");
    Legend11.AddEntry(Match, "Passing Matching", "LPE");
    Legend11.AddEntry(Offline, "Passing Reco", "LPE");
    Legend11.AddEntry(L1_Offline, "Passing L1+Matching+Reco", "LPE");
    Legend11.Draw();

    TLatex Tex1;
    Tex1.SetTextSize(0.03);
    Tex1.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex1.Draw("same");

    TLatex Tex2;
    Tex2.SetTextSize(0.035);
    Tex2.SetTextAlign(31);
    Tex2.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex2.Draw("same");

    c.Write();
    out_file->Close();

    c.SaveAs(output+"/"+DistributionName+".png");
    c.SaveAs(output+"/"+DistributionName+".pdf");
    // c.SetLogy();
    // c.SaveAs(output+"/MuonDistribution_logY.png");
    // c.SaveAs(output+"/MuonDistribution_logY.pdf");
    c.Close();

  }

void PlotTurnOns (TString output, TString TurnOnName, TH1F* L1_Offline, TH1F* Offline, TString YaxisTitle, TString XaxisTitle, Float_t xmin, Float_t xmax, TString type)
  {
    TString rootname = output+"/"+TurnOnName+".root";
    TFile* out_file = TFile::Open(rootname, "recreate");
    out_file->cd();

    TGraphAsymmErrors* TurnOn_L1 = new TGraphAsymmErrors(L1_Offline, Offline, "cp");

    TCanvas c ("c", "c", 900., 900.);
    TPad* p1 = new TPad("p1", "", 0, 0, 1, 1);
    p1->SetGrid();
    p1->SetFillStyle(4000);
    p1->SetFillColorAlpha(0, 0);
    TPad* p2 = new TPad("p2", "", 0, 0, 1, 1);
    p2->SetFillStyle(4000); // will be transparent
    p2->SetFillColorAlpha(0, 0);

    p1->Draw();
    p1->cd();
    TurnOn_L1->SetTitle("");
    TurnOn_L1->GetYaxis()->SetTitle("L1 Efficiency");
    TurnOn_L1->GetXaxis()->SetLimits(xmin,xmax);
    TurnOn_L1->GetXaxis()->SetRangeUser(xmin,xmax);
    TurnOn_L1->GetYaxis()->SetRangeUser(0,1.2);
    TurnOn_L1->GetXaxis()->SetLabelColor(kWhite);
    TurnOn_L1->Draw("ALP");
    TurnOn_L1->SetMarkerStyle(8);
    TurnOn_L1->SetMarkerSize(0.75);
    TurnOn_L1->SetMarkerColor(1);
    TurnOn_L1->SetLineColor(1);
    TurnOn_L1->SetFillColor(3);
    TurnOn_L1->GetHistogram()->SetStats(0);

    Style_t tfont = TurnOn_L1->GetHistogram()->GetYaxis()->GetTitleFont();
    Float_t tsize = TurnOn_L1->GetHistogram()->GetYaxis()->GetTitleSize();
    Style_t lfont = TurnOn_L1->GetHistogram()->GetYaxis()->GetLabelFont();
    Float_t lsize = TurnOn_L1->GetHistogram()->GetYaxis()->GetLabelSize();

    Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
    Double_t ymin = TurnOn_L1->GetMinimum();
    Double_t ymax = TurnOn_L1->GetMaximum();
    Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
    p1->Range(xmin, ymin-0.1*dy, xmax, ymax+0.1*dy);
    p2->Range(xmin, ymin-0.1*dy, xmax, ymax+0.1*dy);
    p2->Draw();
    p2->cd();
    L1_Offline->SetTitle("");
    L1_Offline->Draw("Y+");
    L1_Offline->GetXaxis()->SetTitle(XaxisTitle);
    L1_Offline->GetXaxis()->SetTitleOffset(1.2);
    Int_t steps = (xmax-xmin)/L1_Offline->GetNbinsX();
    TString title = "Entries/"+to_string(steps)+" GeV";
    L1_Offline->GetYaxis()->SetTitle(title);
    L1_Offline->GetXaxis()->SetRangeUser(xmin,xmax);
    L1_Offline->GetYaxis()->SetRangeUser(0, 1.5*L1_Offline->GetMaximum());
    // L1_Offline->GetXaxis()->SetLabelColor(kRed);
    L1_Offline->SetLineWidth(3);
    L1_Offline->SetLineColor(kAzure+5);
    L1_Offline->SetFillColorAlpha(kAzure+5,0.1);
    L1_Offline->SetStats(0);
    Offline->Draw("same");
    Offline->GetYaxis()->SetRangeUser(0, 1.5*L1_Offline->GetMaximum());
    Offline->SetLineWidth(3);
    Offline->SetLineStyle(7);
    Offline->SetLineColor(kAzure);
    Offline->SetStats(0);

    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 50510, "+L");
    axis->SetTitleFont(tfont);
    axis->SetTitleSize(tsize);
    axis->SetLabelFont(lfont);
    axis->SetLabelSize(lsize);
    axis->SetLabelColor(kAzure);
    axis->SetTitleColor(kAzure);
    axis->Draw();
    
    TLatex Tex1;
    Tex1.SetTextSize(0.03);
    Tex1.DrawLatexNDC(0.1,0.91,"#scale[1.5]{CMS} Simulation");
    Tex1.Draw("same");

    TLatex Tex2;
    Tex2.SetTextSize(0.035);
    Tex2.SetTextAlign(31);
    Tex2.DrawLatexNDC(0.9,0.91,"(14 TeV)");
    Tex2.Draw("same");

    TString TurnOn_L1_title;
    TString L1_Offline_title;
    TString Offline_title;

    if (type == "matching")
      {
        TurnOn_L1_title = "Matching Efficiency";
        L1_Offline_title = "Passing Matching";
        Offline_title = "Inclusive";
        L1_Offline->GetYaxis()->SetTitleOffset(1.55);
      }
    
    else if (type == "l1")
      {
        TurnOn_L1_title = "L1 Efficiency";
        L1_Offline_title = "Passing Reco & L1";
        Offline_title = "Passing Reco";
        L1_Offline->GetYaxis()->SetTitleOffset(1.4);
      }

    TLegend *leg = new TLegend(0.55, 0.79, 0.89, 0.89);
    leg->SetFillColor(0);
    leg->SetTextFont(lfont);
    leg->SetTextSize(0.025);
    // leg->SetTextAlign(22);
    leg->AddEntry(TurnOn_L1, TurnOn_L1_title, "L");
    leg->AddEntry(Offline, Offline_title, "L");
    leg->AddEntry(L1_Offline, L1_Offline_title, "L");
    leg->Draw();

    TurnOn_L1->Write();
    Offline->Write();
    L1_Offline->Write();
    c.Write();

    out_file->Close();

    c.SaveAs(output+"/"+TurnOnName+".png");
    c.SaveAs(output+"/"+TurnOnName+".png");
    c.Close();
  }

// Plot jetPt and mjj distribution for VBF events
void PlotEfficiency (TString EventSample)
  {

    TString Method = "Way1"; // when not specified it's way1
    TString JetIDType = "tightLepVetoJetID"; // when not specified it's tightLepVetoJetID
    TString JetSel30 = "JetSel30"; // when not specified it's JetSel30
    bool checkJetSel30 = false;
    if (JetSel30 == "JetSel30") {checkJetSel30 = true;}
    else if (JetSel30 == "NoJetSel30") {checkJetSel30 = false;}

    // TString EventSample = "MC_MiniAOD_EWKino_DemocraticSlepton_mChargino_100to150_17_11_22";
    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/"+EventSample;
    TString output_matching = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotMatchingEfficiency/"+EventSample+"_MatchingEfficiency";
    TString output_L1 = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotMatchingEfficiency/"+EventSample+"_L1Efficiency";

    system("mkdir -p "+output_matching);
    system("mkdir -p "+output_L1);

    // Best L1T selections found in PlotAcceptanceNew/PlotAcceptance_ptjets.C for ptJet1, ptJet2, Mjj, muonPt
    vector<array<Float_t, 4>> set_of_on_cuts;
    vector<array<Float_t, 4>> set_of_off_cuts;

    array<Float_t, 4> trig_on = {70, 30, 220, 6};
    array<Float_t, 4> trig_off = {70+15, 30+15, 220*1.2, 6+2};

    set_of_on_cuts.push_back(trig_on);
    set_of_off_cuts.push_back(trig_off);

    UInt_t nbins_muon = 30;
    UInt_t min_muon = 0;
    UInt_t max_muon = 100;
    TH1F* Total_Muon = new TH1F("Total_Muon", "Total_Muon", nbins_muon, min_muon, max_muon);
    TH1F* Match_Muon = new TH1F("Match_Muon", "Match_Muon", nbins_muon, min_muon, max_muon);
    TH1F* Offline_Muon = new TH1F("Offline_Muon", "Offline_Muon", nbins_muon, min_muon, max_muon);
    TH1F* L1_Offline_Muon = new TH1F("L1_Offline_Muon", "L1_Offline_Muon", nbins_muon, min_muon, max_muon);

    UInt_t nbins_jet1 = 30;
    UInt_t min_jet1 = 0;
    UInt_t max_jet1 = 400;
    TH1F* Total_Jet1 = new TH1F("Total_Jet1", "Total_Jet1", nbins_jet1, min_jet1, max_jet1);
    TH1F* Match_Jet1 = new TH1F("Match_Jet1", "Match_Jet1", nbins_jet1, min_jet1, max_jet1);
    TH1F* Offline_Jet1 = new TH1F("Offline_Jet1", "Offline_Jet1", nbins_jet1, min_jet1, max_jet1);
    TH1F* L1_Offline_Jet1 = new TH1F("L1_Offline_Jet1", "L1_Offline_Jet1", nbins_jet1, min_jet1, max_jet1);

    UInt_t nbins_jet2 = 30;
    UInt_t min_jet2 = 0;
    UInt_t max_jet2 = 200;
    TH1F* Total_Jet2 = new TH1F("Total_Jet2", "Total_Jet2", nbins_jet2, min_jet2, max_jet2);
    TH1F* Match_Jet2 = new TH1F("Match_Jet2", "Match_Jet2", nbins_jet2, min_jet2, max_jet2);
    TH1F* Offline_Jet2 = new TH1F("Offline_Jet2", "Offline_Jet2", nbins_jet2, min_jet2, max_jet2);
    TH1F* L1_Offline_Jet2 = new TH1F("L1_Offline_Jet2", "L1_Offline_Jet2", nbins_jet2, min_jet2, max_jet2);

    UInt_t nbins_mjj = 30;
    UInt_t min_mjj = 0;
    UInt_t max_mjj = 3000;
    TH1F* Total_Mjj = new TH1F("Total_Mjj", "Total_Mjj", nbins_mjj, min_mjj, max_mjj);
    TH1F* Match_Mjj = new TH1F("Match_Mjj", "Match_Mjj", nbins_mjj, min_mjj, max_mjj);
    TH1F* Offline_Mjj = new TH1F("Offline_Mjj", "Offline_Mjj", nbins_mjj, min_mjj, max_mjj);
    TH1F* L1_Offline_Mjj = new TH1F("L1_Offline_Mjj", "L1_Offline_Mjj", nbins_mjj, min_mjj, max_mjj);

    UInt_t events_Total = 0;

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
        // for (UInt_t i_ev = 0 ; i_ev < 10000 ; ++i_ev)
          {
            CheckVBF(inTree, i_ev, set_of_on_cuts, set_of_off_cuts, Total_Muon, Match_Muon, Offline_Muon, L1_Offline_Muon, Total_Jet1, Match_Jet1, Offline_Jet1, L1_Offline_Jet1, Total_Jet2, Match_Jet2, Offline_Jet2, L1_Offline_Jet2, Total_Mjj, Match_Mjj, Offline_Mjj, L1_Offline_Mjj, JetSel30, JetIDType, Method);
          }
      }

    // Plot turn on

    PlotTurnOns (output_matching, "TurnOn_Matching_Muon", Match_Muon, Match_Muon, "Muon Matching Efficiency", "p_{T} Muon [GeV]", min_muon, max_muon, "matching");
    PlotTurnOns (output_matching, "TurnOn_Matching_Jet1", Match_Jet1, Total_Jet1, "Jet1 Matching Efficiency", "p_{T}^{Jet1} [GeV]", min_jet1, max_jet1, "matching");
    PlotTurnOns (output_matching, "TurnOn_Matching_Jet2", Match_Jet2, Total_Jet2, "Jet2 Matching Efficiency", "p_{T}^{Jet2} [GeV]", min_jet2, max_jet2, "matching");
    PlotTurnOns (output_matching, "TurnOn_Matching_Mjj_all", Match_Mjj, Total_Mjj, "Mjj Matching Efficiency", "m_{jj} [GeV]", min_mjj, max_mjj, "matching");

    PlotTurnOns (output_L1, "TurnOn_L1_Muon_all", L1_Offline_Muon, Offline_Muon, "Muon L1 Efficiency", "p_{T} Muon [GeV]", min_muon, max_muon, "l1");
    PlotTurnOns (output_L1, "TurnOn_L1_Jet1_all", L1_Offline_Jet1, Offline_Jet1, "Jet1 L1 Efficiency", "p_{T}^{Jet1} [GeV]", min_jet1, max_jet1, "l1");
    PlotTurnOns (output_L1, "TurnOn_L1_Jet2_all", L1_Offline_Jet2, Offline_Jet2, "Jet2 L1 Efficiency", "p_{T}^{Jet2} [GeV]", min_jet2, max_jet2, "l1");
    PlotTurnOns (output_L1, "TurnOn_L1_Mjj_all",  L1_Offline_Mjj, Offline_Mjj, "Mjj L1 Efficiency", "m_{jj} [GeV]", min_mjj, max_mjj, "l1");

    // Plot distributions

    PlotDistributions (output_matching, Total_Muon, Match_Muon, Offline_Muon, L1_Offline_Muon, "MuonDistribution", "p_{T}^{Offline} (muon) [GeV]");
    PlotDistributions (output_matching, Total_Jet1, Match_Jet1, Offline_Jet1, L1_Offline_Jet1, "Jet1Distribution", "p_{T}^{Offline}(jet1) [GeV]");
    PlotDistributions (output_matching, Total_Jet2, Match_Jet2, Offline_Jet2, L1_Offline_Jet2, "Jet2Distribution", "p_{T}^{Offline}(jet2) [GeV]");
    PlotDistributions (output_matching, Total_Mjj, Match_Mjj, Offline_Mjj, L1_Offline_Mjj, "MjjDistribution", "m_{jj}^{Offline} [GeV]");
    
  }
