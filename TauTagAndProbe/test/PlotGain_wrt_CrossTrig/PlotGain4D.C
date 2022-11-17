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

// Check if a given event passes selections for the acceptance: offline + online + matching for all objects involved in MuTau trigger
bool CheckMuTau (TTree* inTree, UInt_t i_ev)
  {

    ULong64_t       in_EventNumber =  0;
    Int_t           in_RunNumber =  0;
    Int_t           in_lumi =  0;
    vector<float>   *in_PxJet = 0;
    vector<float>   *in_PyJet = 0;
    vector<float>   *in_PzJet = 0;
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

    bool check_MuTau = false;
    bool check_mu_matching = false;
    bool check_tau_matching = false;

    TLorentzVector myGoodOfflineMuon;
    TLorentzVector myGoodOfflineTau;
    TLorentzVector myGoodOfflineJet1;
    TLorentzVector myGoodOfflineJet2;

    TLorentzVector myGoodOnlineMuon;
    TLorentzVector myGoodOnlineTau;

    // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
    bool offline_objects_MuTau = in_PxJet->size() > 1 && in_MuPt->size() > 0 && in_tauPt->size() > 0 ;
    bool L1_objects_MuTau = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 && in_l1tauPt->size() > 0 ;

    if (offline_objects_MuTau && L1_objects_MuTau)
      {
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

        TLorentzVector myGoodOfflineJet1;
        TLorentzVector myGoodOfflineJet2;

        Float_t OfflineJet1Energy = sqrt(pow(in_PxJet->at(0),2)+pow(in_PyJet->at(0),2)+pow(in_PzJet->at(0),2));
        myGoodOfflineJet1.SetPxPyPzE(in_PxJet->at(0), in_PyJet->at(0), in_PzJet->at(0), OfflineJet1Energy);
        Float_t OfflineJet2Energy = sqrt(pow(in_PxJet->at(1),2)+pow(in_PyJet->at(1),2)+pow(in_PzJet->at(1),2));
        myGoodOfflineJet2.SetPxPyPzE(in_PxJet->at(1), in_PyJet->at(1), in_PzJet->at(1), OfflineJet2Energy);

        if (check_mu_matching && check_tau_matching)
          {
            bool check_mu_online    = myGoodOnlineMuon.Pt() > 18;
            bool check_mu_offline   = myGoodOfflineMuon.Pt() > 20;
            bool check_tau_online   = myGoodOnlineTau.Pt() > 24;
            bool check_tau_offline  = myGoodOfflineTau.Pt() > 28;
            bool check_jet1_offline = myGoodOfflineJet1.Pt() > 30;
            bool check_jet2_offline = myGoodOfflineJet2.Pt() > 30;

            bool check_MuTau_online = check_mu_online && check_tau_online;
            bool checkMuTau_offline = check_mu_offline && check_tau_offline && check_jet1_offline && check_jet2_offline;

            check_MuTau = check_MuTau_online && checkMuTau_offline;
          }
      }

    return check_MuTau;
  }


// Check if a given event passes selections for the acceptance: offline + online + matching for all objects involved in VBF trigger
void CheckVBF (TTree* inTree, UInt_t  i_ev, vector<array<Float_t, 4>> set_of_on_cuts, vector<array<Float_t, 4>> set_of_off_cuts, bool pass_MuTau, vector<UInt_t>* acceptance_VBF, vector<UInt_t>* acceptance_MuTau_VBF)
  {

    ULong64_t       in_EventNumber =  0;
    Int_t           in_RunNumber =  0;
    Int_t           in_lumi =  0;
    vector<float>   *in_PxJet = 0;
    vector<float>   *in_PyJet = 0;
    vector<float>   *in_PzJet = 0;
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

    // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
    bool offline_objects_VBF = in_PxJet->size() > 1 && in_MuPt->size() > 0 && in_tauPt->size() > 0 ;
    bool L1_objects_VBF = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 ;

    if (offline_objects_VBF && L1_objects_VBF)
      {

        for (UInt_t i_cut = 0; i_cut < set_of_on_cuts.size(); ++i_cut)
          {

            bool check_VBF = false;
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

            for (UInt_t i_jet1 = 0 ; i_jet1 < in_PxJet->size() ; ++i_jet1)
              {
                Float_t jet1_energy = sqrt(pow(in_PxJet->at(i_jet1),2)+pow(in_PyJet->at(i_jet1),2)+pow(in_PzJet->at(i_jet1),2));
                TLorentzVector myOfflineJet1;
                myOfflineJet1.SetPxPyPzE(in_PxJet->at(i_jet1), in_PyJet->at(i_jet1), in_PzJet->at(i_jet1), jet1_energy);
                for (UInt_t i_jet2 = i_jet1 + 1 ; i_jet2 < in_PxJet->size()-1 ; ++i_jet2)
                  {
                    Float_t jet2_energy = sqrt(pow(in_PxJet->at(i_jet2),2)+pow(in_PyJet->at(i_jet2),2)+pow(in_PzJet->at(i_jet2),2));
                    TLorentzVector myOfflineJet2;
                    myOfflineJet2.SetPxPyPzE(in_PxJet->at(i_jet2), in_PyJet->at(i_jet2), in_PzJet->at(i_jet2), jet2_energy);
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

// Takes 4D rate histogram computed by MakeRatesNew/Rate_ZeroBias_Run316216_new.C and considers only the good combinations giving a rate close to 1 kHz
// The function returns a pointer to the list of online (set_of_on_cuts) and corresponding offline (set_of_off_cuts) cuts giving the good rate
void FindGoodRates(THnF* rate_4D, vector<array<Float_t, 4>>* set_of_on_cuts, vector<array<Float_t, 4>>* set_of_off_cuts, vector<array<Int_t, 4>>* set_of_on_bins, Int_t starting_x_bin)
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

// Loops over the VBF events in MC_MiniAOD_multipleTaus_15_12_21 ntuples and computes the VBF, MuTau and common acceptance
// It returns 3D plots, one for each value of muonPt, of the acceptance of the VBF trigger as a function of the other 3 dimesnions
// It computes the VBF trigger values giving the maximum acceptance and the corresponding gain
// void PlotAcceptance(Int_t starting_x_bin)
void PlotGain()
  {

    Int_t starting_x_bin = 1;

    TString FileName_rate = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_3_1/src/TauTagAndProbe/TauTagAndProbe/test/PureRateStudies/MakePureRates/rate_4D.root";
    TFile f (FileName_rate.Data(),"READ");
    THnF* rate_4D = (THnF*)f.Get("rate_4D");

    vector<array<Float_t, 4>> set_of_on_cuts;
    vector<array<Float_t, 4>> set_of_off_cuts;
    vector<array<Int_t, 4>> set_of_on_bins;

    // definition of online and offline selections list giving a good rate close to 1 kHz 
    FindGoodRates(rate_4D, &set_of_on_cuts, &set_of_off_cuts, &set_of_on_bins, starting_x_bin);

    cout << "Size of the array with Online cuts = " << set_of_on_cuts.size() << endl;
    cout << "Size of the array with Offline cuts = " << set_of_off_cuts.size() << endl;

    UInt_t events_Total = 0; 
    UInt_t acceptance_MuTau = 0;
    vector<UInt_t> acceptance_VBF;
    vector<UInt_t> acceptance_MuTau_VBF;

    cout << "Initialize Acceptance arrays" << endl;
    for (UInt_t i_cut = 0; i_cut < set_of_on_cuts.size(); ++i_cut)
      {
        acceptance_VBF.push_back(0.);
        acceptance_MuTau_VBF.push_back(0.);
      }

    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_3_1/src/TauTagAndProbe/TauTagAndProbe/test/MC_MiniAOD_multipleTaus_05_09_22/";

    // loop over all events in the ntuples to check which one passes the VBF or MuTau selections
    for (int nt = 0 ; nt < 5 ; ++nt)
      {
        TString FileName_ntuples = Path_ntuples + "Ntuple_" + to_string(nt) + ".root";
        TFile f_in(FileName_ntuples.Data(), "READ");
        TDirectoryFile* df = (TDirectoryFile*)f_in.Get("Ntuplizer_noTagAndProbe_multipleTaus");
        TTree* inTree = (TTree*)df->Get("TagAndProbe");

        cout << "\nNumber of events in the Ntuple " << nt << " = " << inTree->GetEntries() << endl;
        events_Total += inTree->GetEntries();

        // each event is analysed for all the possible good combinations of selections
        for (UInt_t i_ev = 0 ; i_ev < inTree->GetEntries() ; ++i_ev)
        // for (UInt_t i_ev = 0 ; i_ev < 20 ; ++i_ev)
          {
            // The MuTau trigger is computed separately, since it has fixed thresholds values
            bool pass_MuTau = CheckMuTau(inTree, i_ev);
            if (pass_MuTau)
              {
                ++ acceptance_MuTau;
              }
            // The VBF trigger is computed for all the possible good combinations of selections
            CheckVBF(inTree, i_ev, set_of_on_cuts, set_of_off_cuts, pass_MuTau, &acceptance_VBF, &acceptance_MuTau_VBF);

          }
      }

    Int_t bins[4] = {rate_4D->GetAxis(0)->GetNbins(), rate_4D->GetAxis(1)->GetNbins(), rate_4D->GetAxis(2)->GetNbins(), rate_4D->GetAxis(3)->GetNbins()};
    Double_t xmin[4] = {rate_4D->GetAxis(0)->GetBinLowEdge(1), rate_4D->GetAxis(1)->GetBinLowEdge(1), rate_4D->GetAxis(2)->GetBinLowEdge(1), rate_4D->GetAxis(3)->GetBinLowEdge(1)};
    Double_t xmax[4] = {rate_4D->GetAxis(0)->GetBinUpEdge(bins[0]), rate_4D->GetAxis(1)->GetBinUpEdge(bins[1]), rate_4D->GetAxis(2)->GetBinUpEdge(bins[2]), rate_4D->GetAxis(3)->GetBinUpEdge(bins[3])};
    // Int_t bins[4] = {20, 20, 32, 13};
    // Double_t xmin[4] = {0., 0., 200., 3.};
    // Double_t xmax[4] = {100., 100., 1000., 16.};

    // inizialization of a vector of histograms, one for each muonPt value
    vector <TH3F*> PtJet1_PtJet2_Mjj;
    for (Int_t m = 0; m < bins[3]; ++m)
      {
        PtJet1_PtJet2_Mjj.push_back(nullptr);
      }

    for (Int_t muon_pt = xmin[3]; muon_pt < xmax[3]; ++muon_pt)
      {
        int m = muon_pt - xmin[3];
        TString histname = Form("PtJet1_PtJet2_Mjj_PtMuon%d", muon_pt);
        PtJet1_PtJet2_Mjj.at(m) = new TH3F(histname, histname, bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1], bins[2], xmin[2], xmax[2]);
      }

    // each histogram is filled with the VBF acceptance corresponding to each combination of jetPt1, jetPt2 and mjj cuts
    for (UInt_t i_c = 0; i_c < set_of_on_cuts.size(); ++i_c)
      {
        for (Int_t muon_pt = xmin[3]; muon_pt < xmax[3]; ++muon_pt)
          {
            if (set_of_on_cuts.at(i_c)[3] == muon_pt)
              {
                int m = muon_pt - xmin[3];
                float gain = Float_t(acceptance_VBF.at(i_c) - acceptance_MuTau_VBF.at(i_c))/Float_t(acceptance_MuTau);
                PtJet1_PtJet2_Mjj.at(m)->SetBinContent(set_of_on_bins.at(i_c)[0], set_of_on_bins.at(i_c)[1], set_of_on_bins.at(i_c)[2], gain);
              }
          }
      }

    // gStyle->SetPalette(kCool);
    TColor::InvertPalette();

    // plot of the histograms
    for (Int_t muon_pt = xmin[3]; muon_pt < xmax[3]; ++muon_pt)
      {
        int m = muon_pt - xmin[3];
        TCanvas c("c","c",900.,650.);
        gPad->SetPad(0.0,0.0,1.0,1.0);
        gPad->SetRightMargin(0.11);
        PtJet1_PtJet2_Mjj.at(m)->SetStats(0);
        PtJet1_PtJet2_Mjj.at(m)->SetTitle("");
        PtJet1_PtJet2_Mjj.at(m)->Draw("BOX2 Z");
        PtJet1_PtJet2_Mjj.at(m)->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
        PtJet1_PtJet2_Mjj.at(m)->GetXaxis()->SetTitleOffset(1.8);
        PtJet1_PtJet2_Mjj.at(m)->GetYaxis()->SetTitle("p_{T}^{j2} [GeV]");
        PtJet1_PtJet2_Mjj.at(m)->GetYaxis()->SetTitleOffset(1.8);
        PtJet1_PtJet2_Mjj.at(m)->GetZaxis()->SetTitle("M_{jj} [GeV]");
        PtJet1_PtJet2_Mjj.at(m)->GetZaxis()->SetTitleOffset(1.4);
        PtJet1_PtJet2_Mjj.at(m)->SetMinimum(0);
        PtJet1_PtJet2_Mjj.at(m)->SetMaximum(0.2);

        TLatex Tex1;
        Tex1.SetTextSize(0.03);
        Tex1.DrawLatexNDC(0.04,0.96,"#scale[1.5]{CMS} Simulation");
        Tex1.Draw("same");

        TLatex Tex2;
        Tex2.SetTextSize(0.035);
        Tex2.SetTextAlign(31);
        Tex2.DrawLatexNDC(0.96,0.96,"(14 TeV)");
        Tex2.Draw("same");

        TLatex Tex3;
        Tex3.SetTextSize(0.035);
        Tex3.SetTextAlign(21);
        Tex3.DrawLatexNDC(0.5,0.93,Form("Gain (p_{T}^{#mu} > %d [GeV])", muon_pt));
        Tex3.Draw("same");

        gPad->Update();
        TPaletteAxis* pal = (TPaletteAxis*) PtJet1_PtJet2_Mjj.at(m)->GetListOfFunctions()->FindObject("palette");
        pal->SetLabelSize(0.01);
        pal->SetX1NDC(0.916); // It doesn't work (Why??)
        PtJet1_PtJet2_Mjj.at(m)->SaveAs(Form("Gain_PtJet1_PtJet2_Mjj_PtMuon%i.root", muon_pt));
        c.SaveAs(Form("Gain_PtJet1_PtJet2_Mjj_PtMuon%d.png", muon_pt));
        c.SaveAs(Form("Gain_PtJet1_PtJet2_Mjj_PtMuon%d.pdf", muon_pt));
      }

    // find the maximum VBF acceptance value and the corresponding selections and gain
    // gain is computed as the number of the events passing only the VBF trigger (and not the MuTau trigger),
    //  with respect to the number of events passing only the MuTau trigger
    UInt_t Max_Acceptance = 0.;
    Float_t Max_Acceptance_gain = 0.;
    array<Float_t, 4> Max_Acceptance_cut;

    for (UInt_t i_c = 0; i_c < set_of_on_cuts.size(); ++i_c)
      {
        Float_t Acceptance_gain = Float_t(acceptance_VBF.at(i_c) - acceptance_MuTau_VBF.at(i_c))/Float_t(acceptance_MuTau);
        if (Acceptance_gain > 0.11 && set_of_on_cuts.at(i_c)[2] < 500)
          {
            cout << "Acceptance cut [" << set_of_on_cuts.at(i_c)[0] << "," << set_of_on_cuts.at(i_c)[1] << "," << set_of_on_cuts.at(i_c)[2] << "," << set_of_on_cuts.at(i_c)[3] << "] VBF = " << acceptance_VBF.at(i_c) << " VBF+MuTau = " << acceptance_MuTau_VBF.at(i_c) << endl;
            cout << "Gain = " << Float_t(acceptance_VBF.at(i_c) - acceptance_MuTau_VBF.at(i_c))/Float_t(acceptance_MuTau) << endl;
          }
        if (Acceptance_gain > Max_Acceptance_gain)
          {
            Max_Acceptance_gain = Acceptance_gain;
            Max_Acceptance_cut = {set_of_on_cuts.at(i_c)[0], set_of_on_cuts.at(i_c)[1], set_of_on_cuts.at(i_c)[2], set_of_on_cuts.at(i_c)[3]};
            Max_Acceptance = acceptance_VBF.at(i_c);
          }
      }

    cout << "Maximum gain is for cut [" << Max_Acceptance_cut[0] << "," << Max_Acceptance_cut[1] << "," << Max_Acceptance_cut[2] << "," << Max_Acceptance_cut[3] << "] = " << Max_Acceptance_gain << endl;
    cout << "Corresponding acceptance = " << Max_Acceptance << endl;

// Plot gain 2D for fixed ptj2 and ptmu

  TH2F* gain_2D = new TH2F("gain_2D","gain_2D", bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2]);

  for (UInt_t i_c = 0; i_c < set_of_on_cuts.size(); ++i_c)
    {
      if (set_of_on_cuts.at(i_c)[1] == 30 && set_of_on_cuts.at(i_c)[3] == 6)
        {
          float gain = Float_t(acceptance_VBF.at(i_c) - acceptance_MuTau_VBF.at(i_c))/Float_t(acceptance_MuTau);
          gain_2D->SetBinContent(set_of_on_bins.at(i_c)[0], set_of_on_bins.at(i_c)[2], gain);
        }
    }

  TCanvas* c_2D = new TCanvas("c_2D","c_2D",700.,550.);
  c_2D->cd();
  gPad->SetPad(0.0,0.0,1.0,1.0);
  gPad->SetTopMargin(2);
  gPad->SetBottomMargin(2);
  gPad->SetRightMargin(2);
  gain_2D->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
  gain_2D->GetXaxis()->SetTitleOffset(1.3);
  gain_2D->GetYaxis()->SetTitle("M_{jj} [GeV]");
  gain_2D->GetYaxis()->SetTitleOffset(1.3);
  gain_2D->GetZaxis()->SetTitle("gain");
  gain_2D->GetZaxis()->SetTitleOffset(0.6);
  gain_2D->SetTitle("");
  gain_2D->SetStats(0);
  gain_2D->Draw("COLZ");
  gain_2D->SetMinimum(0);
  gain_2D->SetMaximum(0.25);

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
  Tex33.DrawLatexNDC(0.25,0.96,"Gain p_{T}^{jet2} > 30 GeV && p_{T}^{#mu} > 6 GeV");
  Tex33.Draw("same");

  c_2D->SaveAs("gain_ptj1_X_ptj2_30_mjj_Y_ptmu_6.png");
  c_2D->SaveAs("gain_ptj1_X_ptj2_30_mjj_Y_ptmu_6.pdf");
  c_2D->Close();


  }