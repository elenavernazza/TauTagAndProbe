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
#include "CheckL1Triggers.h"

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

// void ReadOneVBFEvents (TTree* inTree, UInt_t  i_ev, TString JetIDType, TString Method, bool JetSel30, TLorentzVector* myGoodOfflineTau,
//                   TLorentzVector* myGoodOfflineMuon, TLorentzVector* myGoodOfflineJet1, TLorentzVector* myGoodOfflineJet2, TLorentzVector* myGoodOfflineDiJet, 
//                   TLorentzVector* myGoodOnlineMuon, TLorentzVector* myGoodOnlineJet1, TLorentzVector* myGoodOnlineJet2, TLorentzVector* myGoodOnlineDiJet)
//   {

//     ULong64_t       in_EventNumber =  0;
//     Int_t           in_RunNumber =  0;
//     Int_t           in_lumi =  0;
//     vector<float>   *in_PxJet = 0;
//     vector<float>   *in_PyJet = 0;
//     vector<float>   *in_PzJet = 0;
//     vector<float>   *in_IDJet = 0;
//     vector<float>   *in_l1tPtJet = 0;
//     vector<float>   *in_l1tEtaJet = 0;
//     vector<float>   *in_l1tPhiJet = 0;
//     vector<float>   *in_MuPt = 0;
//     vector<float>   *in_MuEta = 0;
//     vector<float>   *in_MuPhi = 0;
//     vector<bool>    *in_MuIso = 0;
//     vector<int>     *in_MuID = 0;
//     vector<float>   *in_l1tMuPt = 0;
//     vector<float>   *in_l1tMuEta = 0;
//     vector<float>   *in_l1tMuPhi = 0;
//     vector<float>   *in_l1tMuQual = 0;
//     vector<float>   *in_tauPt = 0;
//     vector<float>   *in_tauEta = 0;
//     vector<float>   *in_tauPhi = 0;
//     vector<float>   *in_l1tauPt = 0;
//     vector<float>   *in_l1tauEta = 0;
//     vector<float>   *in_l1tauPhi = 0;

//     inTree->SetBranchAddress("EventNumber", &in_EventNumber);
//     inTree->SetBranchAddress("RunNumber", &in_RunNumber);
//     inTree->SetBranchAddress("lumi", &in_lumi);
//     inTree->SetBranchAddress("jets_px", &in_PxJet);
//     inTree->SetBranchAddress("jets_py", &in_PyJet);
//     inTree->SetBranchAddress("jets_pz", &in_PzJet);
//     inTree->SetBranchAddress("jets_ID", &in_IDJet);
//     inTree->SetBranchAddress("l1tPtJet", &in_l1tPtJet);
//     inTree->SetBranchAddress("l1tEtaJet", &in_l1tEtaJet);
//     inTree->SetBranchAddress("l1tPhiJet", &in_l1tPhiJet);
//     inTree->SetBranchAddress("muons_pt", &in_MuPt);
//     inTree->SetBranchAddress("muons_eta", &in_MuEta);
//     inTree->SetBranchAddress("muons_phi", &in_MuPhi);
//     inTree->SetBranchAddress("muons_PFIsoTight", &in_MuIso);
//     inTree->SetBranchAddress("muons_type", &in_MuID);
//     inTree->SetBranchAddress("l1t_muons_pt", &in_l1tMuPt);
//     inTree->SetBranchAddress("l1t_muons_eta", &in_l1tMuEta);
//     inTree->SetBranchAddress("l1t_muons_phi", &in_l1tMuPhi);
//     inTree->SetBranchAddress("l1t_muons_qual", &in_l1tMuQual);
//     inTree->SetBranchAddress("tauPt", &in_tauPt);
//     inTree->SetBranchAddress("tauEta", &in_tauEta);
//     inTree->SetBranchAddress("tauPhi", &in_tauPhi);
//     inTree->SetBranchAddress("l1tPt", &in_l1tauPt);
//     inTree->SetBranchAddress("l1tEta", &in_l1tauEta);
//     inTree->SetBranchAddress("l1tPhi", &in_l1tauPhi);

//     inTree->GetEntry(i_ev);

//     // remove taus that are matched with jets within deltaR < 0.5
//     vector <bool> jet_is_matched_to_tau;
//     jet_is_matched_to_tau.clear();

//     for (UInt_t ijet = 0; ijet < in_PxJet->size(); ++ijet)
//       {
//         TLorentzVector local_jet;
//         Float_t jet_energy = sqrt(pow(in_PxJet->at(ijet),2)+pow(in_PyJet->at(ijet),2)+pow(in_PzJet->at(ijet),2));
//         local_jet.SetPxPyPzE(in_PxJet->at(ijet), in_PyJet->at(ijet), in_PzJet->at(ijet), jet_energy);

//         bool isMatched = false;
//         for (UInt_t itau = 0; itau < in_tauPt->size(); ++itau)
//           {
//             TLorentzVector local_tau;
//             local_tau.SetPtEtaPhiM(in_tauPt->at(itau), in_tauEta->at(itau), in_tauPhi->at(itau), 0.);
//             if (local_tau.Pt()<10.) continue;
//             if (local_tau.DeltaR(local_jet) < 0.5)
//               {
//                 isMatched = true;
//                 break;
//               }
//           }
//         jet_is_matched_to_tau.push_back(isMatched);
//       }

//     // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
//     bool offline_objects_VBF = in_PxJet->size() > 1 && in_MuPt->size() > 0;
//     bool offline_objects_MuTau = in_tauPt->size() > 0 ;
//     bool L1_objects_VBF = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 ;

//     // if (offline_objects_VBF && L1_objects_VBF && offline_objects_MuTau)
//     //   {

//     bool check_VBF = false; // this has to be inside the loop
//     bool check_VBF_online = false;
//     bool check_VBF_offline = false;
//     bool check_mu_matching = false;
//     bool check_jet1_matching = false;
//     bool check_jet2_matching = false;

//     myGoodOfflineTau->SetPtEtaPhiM(in_tauPt->at(0), in_tauEta->at(0), in_tauPhi->at(0), 1.776);

    // for (UInt_t i_mu = 0 ; i_mu < in_MuPt->size() ; ++i_mu)
    //   {
    //     if (check_mu_matching) break;
    //     if (CheckGoodMuon(in_MuIso->at(i_mu), in_MuID->at(i_mu))) // Isolation and identification of muons
    //       {
    //         TLorentzVector myOfflineMuon;
    //         myOfflineMuon.SetPtEtaPhiM(in_MuPt->at(i_mu), in_MuEta->at(i_mu), in_MuPhi->at(i_mu), 0.105);

    //         for (UInt_t i_L1_mu = 0 ; i_L1_mu < in_l1tMuPt->size() ; ++i_L1_mu)
    //           {
    //             if (check_mu_matching) break;
    //             if (CheckMuonQuality(int(in_l1tMuQual->at(i_L1_mu))) == false) continue;
    //             TLorentzVector myOnlineMuon;
    //             myOnlineMuon.SetPtEtaPhiM(in_l1tMuPt->at(i_L1_mu), in_l1tMuEta->at(i_L1_mu), in_l1tMuPhi->at(i_L1_mu), 0.105);
    //             if (myOfflineMuon.DeltaR(myOnlineMuon) < 0.3)
    //               {
    //                 myGoodOfflineMuon = myOfflineMuon;
    //                 myGoodOnlineMuon = myOnlineMuon;
    //                 check_mu_matching = true;
    //                 break;
    //               }
    //           }
    //       }
    //   }

    // Float_t highest_mjj_offline = -99;
    // int myGoodOfflineJet1Index = -1;
    // int myGoodOfflineJet2Index = -1;
    // int myGoodOnlineJet1Index = -1;
    // int myGoodOnlineJet2Index = -1;

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // /////////////////////////////// Way 1 /////////////////////////////////
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    // if (Method == "Way1")
    //   {
    //     // find the two jets giving highest mjj among jets with pt > 30
    //     for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size() ; ++i_jet1)
    //       {
    //         if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
    //         if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;
    //         Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
    //         TLorentzVector myOfflineJet1;
    //         myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
    //         if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}
    //         for (UInt_t i_jet2 = i_jet1 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
    //           {
    //             if (jet_is_matched_to_tau.at(i_jet2) == true) continue;
    //             if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
    //             Float_t jet2_energy = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
    //             TLorentzVector myOfflineJet2;
    //             myOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy);
    //             if (JetSel30) {if (myOfflineJet2.Pt() < 30) continue;}
    //             TLorentzVector myOfflineDiJet;
    //             myOfflineDiJet = myOfflineJet1 + myOfflineJet2;
    //             if (myOfflineDiJet.M() > highest_mjj_offline)
    //               {
    //                 myGoodOfflineJet1 = myOfflineJet1;
    //                 myGoodOfflineJet1Index = i_jet1;
    //                 myGoodOfflineJet2 = myOfflineJet2;
    //                 myGoodOfflineJet2Index = i_jet2;
    //                 myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;
    //                 highest_mjj_offline = myGoodOfflineDiJet.M();
    //               }
    //           }
    //       }

    //     for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
    //       {
    //         if (check_jet1_matching) break;
    //         TLorentzVector myOnlineJet1;
    //         myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);
    //         if (myGoodOfflineJet1.DeltaR(myOnlineJet1) < 0.5)
    //           {
    //             myGoodOnlineJet1 = myOnlineJet1;
    //             myGoodOnlineJet1Index = i_L1_jet1;
    //             check_jet1_matching = true;
    //             break;
    //           }
    //       }

    //     for (UInt_t i_L1_jet2 = 0 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
    //       {
    //         if (check_jet2_matching) break;
    //         if (int(i_L1_jet2) != myGoodOnlineJet1Index)
    //           {
    //             TLorentzVector myOnlineJet2;
    //             myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);
    //             if (myGoodOfflineJet2.DeltaR(myOnlineJet2) < 0.5)
    //               {
    //                 myGoodOnlineJet2 = myOnlineJet2;
    //                 myGoodOnlineJet2Index = i_L1_jet2;
    //                 check_jet2_matching = true;
    //                 break;
    //               }
    //           }
    //       }

    //     myGoodOnlineDiJet = myGoodOnlineJet1 + myGoodOnlineJet2;

    //   }

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // /////////////////////////////// Way 2 /////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // if (Method == "Way2")
    //   {

    //     for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size()-1 ; ++i_jet1)
    //       {
    //         if (check_jet1_matching && check_jet2_matching) break;
    //         if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
    //         if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;

    //         // compute jet 1 energy from px, py, pz information assuming jet mass to 0
    //         Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
    //         // definition of Lorentz vector for offline jet 1 (px, py, pz, E [GeV])
    //         TLorentzVector myOfflineJet1;
    //         myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
    //         if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}

    //         for (UInt_t i_jet2 = i_jet1+1 ; i_jet2 < in_PxJet->size() ; ++i_jet2)
    //           {
    //             if (check_jet1_matching && check_jet2_matching) break;
    //             if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
    //             // compute jet 2 energy from px, py, pz information assuming jet mass to 0
    //             Float_t jet2_energy = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
    //             // definition of Lorentz vector for offline jet 2 (px, py, pz, E [GeV])
    //             TLorentzVector myOfflineJet2;
    //             myOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy);
    //             if (JetSel30) {if (myOfflineJet2.Pt() < 30) continue;}

    //             for (UInt_t i_L1_jet1 = 0 ; i_L1_jet1 < in_l1tPtJet->size() ; ++i_L1_jet1)
    //               {

    //                 if (check_jet1_matching && check_jet2_matching) break;
    //                 // definition of Lorentz vector for online jet 1 (pt, eta, phi, mass [GeV])
    //                 TLorentzVector myOnlineJet1;
    //                 myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet1), in_l1tEtaJet->at(i_L1_jet1), in_l1tPhiJet->at(i_L1_jet1), 0);

    //                 for (UInt_t i_L1_jet2 = i_L1_jet1+1 ; i_L1_jet2 < in_l1tPtJet->size() ; ++i_L1_jet2)
    //                   {

    //                     if (check_jet1_matching && check_jet2_matching) break;
    //                     // definition of Lorentz vector for online jet 2 (pt, eta, phi, mass [GeV])
    //                     TLorentzVector myOnlineJet2;
    //                     myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(i_L1_jet2), in_l1tEtaJet->at(i_L1_jet2), in_l1tPhiJet->at(i_L1_jet2), 0);

    //                     // check matching based on deltaR between online and offline object
    //                     if (myOfflineJet1.DeltaR(myOnlineJet1) < 0.5 && myOfflineJet2.DeltaR(myOnlineJet2) < 0.5)
    //                       {
    //                         myGoodOfflineJet1 = myOfflineJet1;
    //                         myGoodOfflineJet1Index = i_jet1;
    //                         myGoodOfflineJet2 = myOfflineJet2;
    //                         myGoodOfflineJet2Index = i_jet2;
    //                         myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;

    //                         myGoodOnlineJet1 = myOnlineJet1;
    //                         myGoodOnlineJet1Index = i_L1_jet1;
    //                         myGoodOnlineJet2 = myOnlineJet2;
    //                         myGoodOnlineJet2Index = i_L1_jet2;
    //                         myGoodOnlineDiJet = myGoodOnlineJet1 + myGoodOnlineJet2;

    //                         check_jet1_matching = true;
    //                         check_jet2_matching = true;
    //                         break;
    //                       }
    //                   }
    //               }
    //           }
    //       }
    //   }

  // }

// Check if a given event passes selections for the acceptance: offline + online + matching for all objects involved in VBF trigger
void CheckVBF_vs_MuTau (TTree* inTree, UInt_t  i_ev, vector<array<Float_t, 4>> set_of_on_cuts, vector<array<Float_t, 4>> set_of_off_cuts, bool pass_MuTau, vector<UInt_t>* acceptance_VBF, vector<UInt_t>* acceptance_MuTau_VBF, TString JetIDType, TString Method, bool JetSel30)
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
                if (CheckGoodMuon(in_MuIso->at(i_mu), in_MuID->at(i_mu))) // Isolation and identification of muons
                  {
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
                    if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
                    if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;
                    Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                    TLorentzVector myOfflineJet1;
                    myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
                    if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}
                    for (UInt_t i_jet2 = i_jet1 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
                      {
                        if (jet_is_matched_to_tau.at(i_jet2) == true) continue;
                        if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
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
                    if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
                    if (CheckGoodJet(in_IDJet->at(i_jet1), JetIDType) == false) continue;

                    // compute jet 1 energy from px, py, pz information assuming jet mass to 0
                    Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                    // definition of Lorentz vector for offline jet 1 (px, py, pz, E [GeV])
                    TLorentzVector myOfflineJet1;
                    myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
                    if (JetSel30) {if (myOfflineJet1.Pt() < 30) continue;}

                    for (UInt_t i_jet2 = i_jet1+1 ; i_jet2 < in_PxJet->size() ; ++i_jet2)
                      {
                        if (check_jet1_matching && check_jet2_matching) break;
                        if (CheckGoodJet(in_IDJet->at(i_jet2), JetIDType) == false) continue;
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
                bool check_jet1_online = myGoodOnlineJet1.Pt() > set_of_on_cuts.at(i_cut)[0];
                bool check_jet2_online = myGoodOnlineJet2.Pt() > set_of_on_cuts.at(i_cut)[1];
                bool check_mjj_online  = myGoodOnlineDiJet.M() > set_of_on_cuts.at(i_cut)[2]; // this is the problem
                bool check_mu_online   = myGoodOnlineMuon.Pt() > set_of_on_cuts.at(i_cut)[3];
                
                bool check_jet1_offline = myGoodOfflineJet1.Pt() > set_of_off_cuts.at(i_cut)[0];
                bool check_jet2_offline = myGoodOfflineJet2.Pt() > set_of_off_cuts.at(i_cut)[1];
                bool check_mjj_offline  = myGoodOfflineDiJet.M() > set_of_off_cuts.at(i_cut)[2];
                bool check_mu_offline   = myGoodOfflineMuon.Pt() > set_of_off_cuts.at(i_cut)[3];
                bool check_tau_offline   = myGoodOfflineTau.Pt() > 20;

                check_VBF_online = check_jet1_online && check_jet2_online && check_mjj_online && check_mu_online;
                check_VBF_offline = check_jet1_offline && check_jet2_offline && check_mjj_offline && check_mu_offline && check_tau_offline;

                check_VBF = check_VBF_online && check_VBF_offline;

                // if (check_VBF_online)
                //   {
                //     cout << "\nEvent number " << i_ev << endl;
                //     cout << "\nOffline Muon Pt = ";
                //     for (UInt_t a = 0; a < in_MuPt->size(); ++a) cout << "[" << in_MuPt->at(a) << ", " << in_MuEta->at(a) << ", " << in_MuPhi->at(a) << "] ";
                //     cout << "\nOffline Good Muon Pt = ";
                //     for (UInt_t b = 0; b < in_MuPt->size(); ++b) 
                //       {
                //         if (CheckGoodMuon(in_MuIso->at(b), in_MuID->at(b)))
                //           {
                //             cout << "[" << in_MuPt->at(b) << ", " << in_MuEta->at(b) << ", " << in_MuPhi->at(b) << "] ";
                //           }
                //       }
                //   }
              }

            // if the event passes the VBF trigger selections, it is considered in the computation of the acceptance
            if (check_VBF)
              {
                acceptance_VBF->at(i_cut) += 1.;
                if (pass_MuTau)
                  {
                    acceptance_MuTau_VBF->at(i_cut) += 1.;
                  }
              }
          }
      }
  }

// Check if a given event passes selections for the acceptance: offline + online + matching for all objects involved in MuTau trigger
bool CheckMuTau (TTree* inTree, UInt_t i_ev, TString JetIDType, TString Method, bool JetSel30)
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

    bool check_MuTau = false;

    // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
    bool offline_objects_MuTau = in_MuPt->size() > 0 && in_tauPt->size() > 0 ;
    bool offline_objects_VBF = in_PxJet->size() > 1;
    bool L1_objects_MuTau = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 && in_l1tauPt->size() > 0 ;

    if (offline_objects_MuTau && L1_objects_MuTau && offline_objects_VBF)
      {

        TLorentzVector myGoodOfflineMuon;
        TLorentzVector myGoodOfflineTau;
        TLorentzVector myGoodOfflineJet1;
        TLorentzVector myGoodOfflineJet2;
        TLorentzVector myGoodOfflineDiJet;

        TLorentzVector myGoodOnlineMuon;
        TLorentzVector myGoodOnlineTau;

        // Muon definition

        bool check_mu_matching = false;
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
                    if (CheckMuonQuality(int(in_l1tMuQual->at(i_L1_mu))) == false) continue;
                    TLorentzVector myOnlineMuon;
                    myOnlineMuon.SetPtEtaPhiM(in_l1tMuPt->at(i_L1_mu), in_l1tMuEta->at(i_L1_mu), in_l1tMuPhi->at(i_L1_mu), 0.105);
                    Float_t deltaR = myOfflineMuon.DeltaR(myOnlineMuon);
                    if (deltaR < 0.3)
                      {
                        check_mu_matching = true;
                        myGoodOfflineMuon = myOfflineMuon;
                        myGoodOnlineMuon = myOnlineMuon;
                        break;
                      }
                  }
              }
          }

        // Tau definition

        bool check_tau_matching = false;
        for (UInt_t i_tau = 0 ; i_tau < in_tauPt->size() ; ++i_tau)
          {
            if (check_tau_matching) break;
            TLorentzVector myOfflineTau;
            myOfflineTau.SetPtEtaPhiM(in_tauPt->at(i_tau), in_tauEta->at(i_tau), in_tauPhi->at(i_tau), 1.776);

            for (UInt_t i_L1_tau = 0 ; i_L1_tau < in_l1tauPt->size() ; ++i_L1_tau)
              {
                if (check_tau_matching) break;
                TLorentzVector myOnlineTau;
                myOnlineTau.SetPtEtaPhiM(in_l1tauPt->at(i_L1_tau), in_l1tauEta->at(i_L1_tau), in_l1tauPhi->at(i_L1_tau), 1.776);
                Float_t deltaR = myOfflineTau.DeltaR(myOnlineTau);
                if (deltaR < 0.3)
                  {
                    check_tau_matching = true;
                    myGoodOfflineTau = myOfflineTau;
                    myGoodOnlineTau = myOnlineTau;
                    break;
                  }
              }
          }

        // Offline selection on Mjj and jets

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

        if (Method == "Way1")
          {
            Float_t highest_mjj_offline_way1 = -99;
            UInt_t myGoodOfflineJet1Index_way1 = -1;
            UInt_t myGoodOfflineJet2Index_way1 = -1;

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
                        myGoodOfflineJet1 = myOfflineJet1_way1;
                        myGoodOfflineJet2 = myOfflineJet2_way1;
                        myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;
                        highest_mjj_offline_way1 = myGoodOfflineDiJet.M();
                        myGoodOfflineJet1Index_way1 = i_jet1;
                        myGoodOfflineJet2Index_way1 = i_jet2;
                      }
                  }
              }
          }

        if (Method == "Way2")
          {
            // first offline jet not matched to any tau
            int myGoodOfflineJet1Index_way2 = -1;
            for (UInt_t i_jet1 = 0; i_jet1 < in_PxJet->size(); ++i_jet1)
              {
                if (jet_is_matched_to_tau.at(i_jet1) == true) continue;
                  {
                    Float_t OfflineJet1Energy_way2 = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                    myGoodOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), OfflineJet1Energy_way2);
                    break;
                  }
              }
            for (UInt_t i_jet2 = myGoodOfflineJet1Index_way2 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
              {
                if (jet_is_matched_to_tau.at(i_jet2) == true) continue;
                  {
                    Float_t OfflineJet2Energy_way2 = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                    myGoodOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), OfflineJet2Energy_way2);
                    break;            
                  }            
              }
            myGoodOfflineDiJet = myGoodOfflineJet1 + myGoodOfflineJet2;
          }

        // Print jets
        // TLorentzVector myLeadingJet;
        // Float_t myLeadingJetEnergy = sqrt(pow(in_PxJet->at(0),2)+pow(in_PyJet->at(0),2)+pow(in_PzJet->at(0),2));
        // myLeadingJet.SetPxPyPzE(in_PxJet->at(0), in_PyJet->at(0), in_PzJet->at(0), myLeadingJetEnergy);
        // TLorentzVector mySubleadingJet;
        // Float_t mySubleadingJetEnergy = sqrt(pow(in_PxJet->at(1),2)+pow(in_PyJet->at(1),2)+pow(in_PzJet->at(1),2));
        // mySubleadingJet.SetPxPyPzE(in_PxJet->at(1), in_PyJet->at(1), in_PzJet->at(1), mySubleadingJetEnergy);
        // TLorentzVector myFirstDiJet;
        // myFirstDiJet = myLeadingJet + mySubleadingJet;

        // cout << "\nEvent number " << i_ev << endl;
        // cout << "\nOffline Jets Pt = [ ";
        // for (UInt_t a = 0; a < in_PxJet->size(); ++a) 
        //   {
        //     Float_t energy = sqrt(pow(in_PxJet->at(a),2)+pow(in_PyJet->at(a),2)+pow(in_PzJet->at(a),2));
        //     TLorentzVector myJet;
        //     myJet.SetPxPyPzE(in_PxJet->at(a), in_PyJet->at(a), in_PzJet->at(a), energy);

        //     cout << ", " << myJet.Pt() ;
        //   }
        // cout << " ]" << endl;
        // cout << "Mjj (jet1 and jet2) = " << myFirstDiJet.M() << endl;
        // cout << "Mjj (jet" << myGoodOfflineJet1Index << " and jet" << myGoodOfflineJet2Index << ") = " << myGoodOfflineDiJet.M();
        // cout << " DiJet Pt = " << myGoodOfflineDiJet.Pt();
        // cout << " DiJet Pz = " << myGoodOfflineDiJet.Pz() << endl;
        // cout << "Jet" << myGoodOfflineJet1Index << " Pt = " << myGoodOfflineJet1.Pt() << ", Jet" << myGoodOfflineJet2Index << " Pt = " << myGoodOfflineJet2.Pt() << endl;
        // cout << "Jet" << myGoodOfflineJet1Index << " Pz = " << myGoodOfflineJet1.Pz() << ", Jet" << myGoodOfflineJet2Index << " Pz = " << myGoodOfflineJet2.Pz() << endl;

        bool check_mu_online    = myGoodOnlineMuon.Pt() > 18;
        bool check_mu_offline   = myGoodOfflineMuon.Pt() > 20;
        bool check_tau_online   = myGoodOnlineTau.Pt() > 24;
        bool check_tau_offline  = myGoodOfflineTau.Pt() > 28;
        bool check_jet1_offline = myGoodOfflineJet1.Pt() > 30;
        bool check_jet2_offline = myGoodOfflineJet2.Pt() > 30;
        bool check_mjj_offline  = myGoodOfflineDiJet.M() > 200;

        bool check_MuTau_online = check_mu_online && check_tau_online;
        bool check_MuTau_offline = check_mu_offline && check_tau_offline && check_jet1_offline && check_jet2_offline;

        check_MuTau = check_MuTau_online && check_MuTau_offline;
      }

    return check_MuTau;
  }

// Takes 4D rate histogram computed by MakeRatesNew/Rate_ZeroBias_Run316216_new.C and considers only the good combinations giving a rate close to 1 kHz
// The function returns a pointer to the list of online (set_of_on_cuts) and corresponding offline (set_of_off_cuts) cuts giving the good rate
void FindGoodRates (THnF* rate_4D, vector<array<Float_t, 4>>* set_of_on_cuts, vector<array<Float_t, 4>>* set_of_off_cuts, vector<array<Int_t, 4>>* set_of_on_bins, Int_t starting_x_bin)
  {

    cout << "\nFind Good Rates\n" << endl;
    UInt_t combinations = 0;
    Int_t bins[4] = {rate_4D->GetAxis(0)->GetNbins(), rate_4D->GetAxis(1)->GetNbins(), rate_4D->GetAxis(2)->GetNbins(), rate_4D->GetAxis(3)->GetNbins()};

    for (Int_t i_bin_y = starting_x_bin ; i_bin_y <= bins[1] ; ++i_bin_y)
      {
        for (Int_t i_bin_x = i_bin_y ; i_bin_x <= bins[0] ; ++i_bin_x) // To avoid cuts where ptj1 < ptj2
          {
            for (Int_t i_bin_z = 1 ; i_bin_z <= bins[2] ; ++i_bin_z)
              {
                for (Int_t i_bin_w = 1 ; i_bin_w <= bins[3] ; ++i_bin_w)
                  {
                    Int_t v_bins[4] = {i_bin_x, i_bin_y, i_bin_z, i_bin_w};
                    Float_t rate = rate_4D->GetBinContent(v_bins);

                    if (rate_4D->GetBinContent(v_bins) > 0.95 && rate_4D->GetBinContent(v_bins) < 1.05)
                      {
                        Float_t trig_x = Float_t(rate_4D->GetAxis(0)->GetBinLowEdge(i_bin_x));
                        Float_t trig_y = Float_t(rate_4D->GetAxis(1)->GetBinLowEdge(i_bin_y));
                        Float_t trig_z = Float_t(rate_4D->GetAxis(2)->GetBinLowEdge(i_bin_z));
                        Float_t trig_w = Float_t(rate_4D->GetAxis(3)->GetBinLowEdge(i_bin_w));

                        ++ combinations;

                        array<Float_t, 4> trig_on = {trig_x, trig_y, trig_z, trig_w};
                        array<Float_t, 4> trig_off = {Float_t(trig_x+15.), Float_t(trig_y+15.), Float_t(trig_z*1.2), Float_t(trig_w+2.)};
                        array<Int_t, 4> on_bins = {i_bin_x, i_bin_y, i_bin_z, i_bin_w};
                        set_of_on_cuts->push_back(trig_on);
                        set_of_off_cuts->push_back(trig_off);
                        set_of_on_bins->push_back(on_bins);
                      }
                  }
              }
          }
      }
  }

void PlotGain_2D_ptj1_mjj (Float_t fixed_cut_y, Float_t fixed_cut_w, Int_t bin0, Int_t xmin0, Int_t xmax0, Int_t bin2, Int_t xmin2, Int_t xmax2, vector<array<Float_t, 4>> set_of_on_cuts, vector<array<Int_t, 4>> set_of_on_bins, vector<UInt_t> acceptance_VBF, vector<UInt_t> acceptance_MuTau_VBF, UInt_t acceptance_MuTau, TString Output)
  {

    // Plot gain 2D for fixed ptj2 and ptmu
    TH2F* gain_2D_ptj1_mjj = new TH2F("gain_2D_ptj1_mjj","gain_2D_ptj1_mjj", bin0, xmin0, xmax0, bin2, xmin2, xmax2);

    for (UInt_t i_c = 0; i_c < set_of_on_cuts.size(); ++i_c)
      {
        if (set_of_on_cuts.at(i_c)[1] == fixed_cut_y && set_of_on_cuts.at(i_c)[3] == fixed_cut_w)
          {
            float gain = Float_t(acceptance_VBF.at(i_c) - acceptance_MuTau_VBF.at(i_c))/Float_t(acceptance_MuTau);
            gain_2D_ptj1_mjj->SetBinContent(set_of_on_bins.at(i_c)[0], set_of_on_bins.at(i_c)[2], gain);
            // cout << gain << endl;
            // cout << gain_2D_ptj1_mjj->GetBinContent(set_of_on_bins.at(i_c)[0], set_of_on_bins.at(i_c)[2]) << endl;
          }
      }

    gStyle->SetPaintTextFormat(".3f");
    TCanvas* c_2D = new TCanvas("c_2D","c_2D",700.,550.);
    c_2D->cd();
    c_2D->SetRightMargin(0.14);
    gain_2D_ptj1_mjj->GetXaxis()->SetTitle("p_{T}^{j1} > X [GeV]");
    gain_2D_ptj1_mjj->GetXaxis()->SetTitleOffset(1.3);
    gain_2D_ptj1_mjj->GetYaxis()->SetTitle("M_{jj} > Y [GeV]");
    gain_2D_ptj1_mjj->GetYaxis()->SetTitleOffset(1.3);
    gain_2D_ptj1_mjj->GetZaxis()->SetTitle("gain");
    gain_2D_ptj1_mjj->GetZaxis()->SetTitleOffset(1.1);
    gain_2D_ptj1_mjj->SetTitle("");
    gain_2D_ptj1_mjj->SetStats(0);
    gain_2D_ptj1_mjj->Draw("COLZ");
    gain_2D_ptj1_mjj->Draw("text same");
    gain_2D_ptj1_mjj->SetMinimum(0.);
    gain_2D_ptj1_mjj->SetMaximum(max(1.,1.3*gain_2D_ptj1_mjj->GetMaximum()));

    TLatex Tex11;
    Tex11.SetTextSize(0.03);
    Tex11.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex11.Draw("same");

    TLatex Tex22;
    Tex22.SetTextSize(0.035);
    Tex22.SetTextAlign(31);
    Tex22.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex22.Draw("same");

    TLatex Tex33;
    Tex33.SetTextSize(0.035);
    Tex33.SetTextAlign(12);
    TString printout = "Gain p_{T}^{jet2} > " + to_string(int(fixed_cut_y)) + " GeV && p_{T}^{#mu} > " + to_string(int(fixed_cut_w)) + " GeV wrt Mu18Tau24";
    Tex33.DrawLatexNDC(0.25,0.96,printout);
    Tex33.Draw("same");

    c_2D->SaveAs(Output+"/Gain_2D_Mu18Tau24_ptj1_X_ptj2_"+to_string(int(fixed_cut_y))+"_mjj_Y_ptmu_"+to_string(int(fixed_cut_w))+".png");
    c_2D->SaveAs(Output+"/Gain_2D_Mu18Tau24_ptj1_X_ptj2_"+to_string(int(fixed_cut_y))+"_mjj_Y_ptmu_"+to_string(int(fixed_cut_w))+".pdf");
    c_2D->Close();

  }

void PrintMuons (TTree* inTree, UInt_t  i_ev)
  {
    vector<float>   *in_MuPt = 0;
    vector<float>   *in_MuEta = 0;
    vector<float>   *in_MuPhi = 0;
    vector<bool>    *in_MuIso = 0;
    vector<int>     *in_MuID = 0;

    inTree->SetBranchAddress("muons_pt", &in_MuPt);
    inTree->SetBranchAddress("muons_eta", &in_MuEta);
    inTree->SetBranchAddress("muons_phi", &in_MuPhi);
    inTree->SetBranchAddress("muons_PFIsoTight", &in_MuIso);
    inTree->SetBranchAddress("muons_type", &in_MuID);

    inTree->GetEntry(i_ev);

    // cout << "\nEvent number = " << i_ev << "\nOffline Muon Pt = ";
    for (UInt_t a = 0; a < in_MuPt->size(); ++a) 
      {
        if (CheckGoodMuon(in_MuIso->at(a), in_MuID->at(a))) cout << "[" << in_MuPt->at(a) << ", " << in_MuEta->at(a) << ", " << in_MuPhi->at(a) << ", " << in_MuIso->at(a) << ", " << in_MuID->at(a) << "] ";
      }
  }
