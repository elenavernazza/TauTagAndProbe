#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TLorentzVector.h>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>

// for the ZeroBias sample 362616 we can look at this page
// https://cmsoms.cern.ch/cms/runs/lumisection?cms_run=362616&cms_run_sequence=GLOBAL-RUN
// thisLumiRun = 20.5 → 2.050E34;

using namespace std;

void Rate()
{

  Float_t Total = 0;
  Float_t Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = 0.;

  cout << "Begin loop" << endl;

  TString path = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/ZeroBias_Run362616/";

  float nb = 2450.;
  float scale_rate = 0.001*(nb*11245.6);
  float thisLumiRun = 2.050E34;
  float scaleToLumi = 2.000E34;
  float scale_lumi = scaleToLumi/thisLumiRun;

  TString FileName_in = path + "/EphemeralZeroBias0__Run2022G_Run362616__RAW.root" ;
  TFile f_in(FileName_in.Data(),"READ");
  TDirectoryFile* df = (TDirectoryFile*)f_in.Get("ZeroBias");
  TTree* inTree = (TTree*)df->Get("ZeroBias");

  ULong64_t       in_EventNumber =  0;
  Int_t           in_RunNumber =  0;
  Int_t           in_lumi =  0;
  vector<float>   *in_l1tPtJet = 0;
  vector<float>   *in_l1tEtaJet = 0;
  vector<float>   *in_l1tPhiJet = 0;
  vector<vector<float>>   *in_l1tMuPt = 0;
  vector<vector<float>>   *in_l1tMuEta = 0; 
  vector<vector<float>>   *in_l1tMuPhi = 0; 
  vector<vector<int>>     *in_l1tMuQual = 0;

  inTree->SetBranchAddress("EventNumber", &in_EventNumber);
  inTree->SetBranchAddress("RunNumber", &in_RunNumber);
  inTree->SetBranchAddress("lumi", &in_lumi);
  inTree->SetBranchAddress("l1tPtJet", &in_l1tPtJet);
  inTree->SetBranchAddress("l1tEtaJet", &in_l1tEtaJet);
  inTree->SetBranchAddress("l1tPhiJet", &in_l1tPhiJet);
  inTree->SetBranchAddress("l1tMuPt", &in_l1tMuPt);
  inTree->SetBranchAddress("l1tMuEta", &in_l1tMuEta);
  inTree->SetBranchAddress("l1tMuPhi", &in_l1tMuPhi);
  inTree->SetBranchAddress("l1tMuQual", &in_l1tMuQual);

  cout << "\nNumber of events in the Ntuple = " << inTree->GetEntries() << endl;
  for(UInt_t i = 0 ; i < inTree->GetEntries() ; ++i)
    {
      if (i%100000 == 0) cout << "Analysing entry " << i << endl;
      inTree->GetEntry(i);
      scale_lumi = scaleToLumi/thisLumiRun;
      ++Total;

      // at least 2 jets and one muon
      bool L1_objects = in_l1tPtJet->size() > 1 && in_l1tMuPt->at(2).size() > 0 ;

      TLorentzVector myGoodOnlineMuon;
      TLorentzVector myGoodOnlineJet1;
      TLorentzVector myGoodOnlineJet2;
      TLorentzVector myGoodOnlineDiJet;
      float myGoodOnlineMjj = -1.;

      if (L1_objects)
        {

          // to call the muon we have to take the second bunch crossing (in the center)
          UInt_t good_BX = 2;
          // take the first (most enegrgetic) muon passing the open quality requirement
          for (UInt_t iMuon = 0; iMuon < in_l1tMuPt->at(good_BX).size(); ++iMuon)
            {
              if (in_l1tMuQual->at(good_BX).at(iMuon) < 4) continue; // open quality: 4, 5, 6, 7, 8, 9, 10, 11, 12,13, 14, 15
              myGoodOnlineMuon.SetPtEtaPhiM(in_l1tMuPt->at(good_BX).at(iMuon),in_l1tMuEta->at(good_BX).at(iMuon),in_l1tMuPhi->at(good_BX).at(iMuon),0.105);
              break;
            }

          // loop among the pairs of jets (80,30) giving the highest mjj
          for (UInt_t iL1Jet1 = 0 ; iL1Jet1 < in_l1tPtJet->size() ; ++iL1Jet1)
            {
              TLorentzVector myOnlineJet1;
              myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(iL1Jet1),in_l1tEtaJet->at(iL1Jet1),in_l1tPhiJet->at(iL1Jet1),0.);
              if (myOnlineJet1.Pt() < 80) continue;

              for (UInt_t jL1Jet2 = iL1Jet1 + 1 ; jL1Jet2 < in_l1tPtJet->size() ; ++jL1Jet2)
                {
                  TLorentzVector myOnlineJet2;
                  myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(jL1Jet2),in_l1tEtaJet->at(jL1Jet2),in_l1tPhiJet->at(jL1Jet2),0.);
                  if (myOnlineJet2.Pt() < 30) continue;
                  
                  TLorentzVector myOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                  float myOnlineMjj = myOnlineDiJet.M();

                  if (myOnlineMjj > myGoodOnlineMjj)
                    {
                      myGoodOnlineMjj = myOnlineMjj;
                      myGoodOnlineJet1 = myOnlineJet1;
                      myGoodOnlineJet2 = myOnlineJet2;
                      myGoodOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                    }
                }  
            }       

          if (myGoodOnlineJet1.Pt() >= 80 && myGoodOnlineJet2.Pt() >= 30 && myGoodOnlineDiJet.M() >= 420 && myGoodOnlineMuon.Pt() >= 8) 
            {
              Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 += 1.;
            }
        }
    }

  cout << "Total number of events = " << Total << endl;
  cout << "Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 << endl;
  cout << "Rate L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8/Total*scale_rate*scale_lumi*1000 << " Hz" << endl;

  return;

}
