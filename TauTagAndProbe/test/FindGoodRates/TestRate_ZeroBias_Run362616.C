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
#include <cstring>  

// for the ZeroBias sample 362616 we can look at this page
// https://cmsoms.cern.ch/cms/runs/lumisection?cms_run=362616&cms_run_sequence=GLOBAL-RUN
// thisLumiRun = 20.5 → 2.050E34;

using namespace std;

void Rate()
{

  ULong64_t v_events[84] = { 206006627, 177670917, 175861870, 175394993, 232614668, 524787473, 639570227, 700505402, 325211142, 267082889, 218458730, \
                          660053760, 680025275, 345808342, 298800791, 637233481, 724304890, 486436299, 511491518, 674997190, 325185999, 376187696, \
                          692892160, 473916757, 475112445, 330563641, 517757622, 504700767, 208750421, 371351117, 774081902, 411536592, 716908008, \
                          730337285, 466059424, 377602105, 698332779, 335953825, 337574840, 662608795, 538600619, 637644228, 183150718, 331539331, \
                          713391613, 363858525, 593955846, 777970175, 461620678, 594933681, 593430504, 183762179, 222726465, 174996454, 771555726, \
                          361902139, 360580557, 546497022, 646020636, 537261593, 584043008, 294782474, 295577408, 323402859, 382723269, 743556318, \
                          675428368, 213221794, 216253377, 282969787, 320993378, 210603688, 666269100, 271763078, 702922893, 766406432, 456796034, \
                          601219979, 600960634, 575772020, 419108697, 661379920, 363850694, 363690468 };

  cout << v_events[0] << endl;

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
      // if (i%100000 == 0) cout << "Analysing entry " << i << endl;
      inTree->GetEntry(i);
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
              for (int i_fired = 0; i_fired < 84; ++i_fired)
                {
                  if (in_EventNumber == v_events[i_fired])
                    {
                      if (myGoodOnlineMuon.Pt() >= 8) cout << "Found open muon" << endl;
                      break;
                    }
                }
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

          for (int i_fired = 0; i_fired < 84; ++i_fired)
            {
              if (in_EventNumber == v_events[i_fired])
                {
                  cout << in_EventNumber << " " << myGoodOnlineJet1.Pt() << " " << myGoodOnlineJet2.Pt() << " ";
                  cout << myGoodOnlineDiJet.M() << " " << myGoodOnlineMuon.Pt();
                  cout << " " << (myGoodOnlineJet1.Pt() >= 80) << " " << (myGoodOnlineJet2.Pt() >= 30) << " " << (myGoodOnlineDiJet.M() >= 420) << " " << (myGoodOnlineMuon.Pt() >= 8);
                  break;
                }
            }

          if (myGoodOnlineJet1.Pt() >= 80 && myGoodOnlineJet2.Pt() >= 30 && myGoodOnlineDiJet.M() >= 420 && myGoodOnlineMuon.Pt() >= 8) 
            {
              Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 += 1.;
              cout << " passed" << endl;
              // cout << in_EventNumber << " " << myGoodOnlineJet1.Pt() << " " << myGoodOnlineJet2.Pt() << " ";
              // cout << myGoodOnlineDiJet.M() << " " << myGoodOnlineMuon.Pt() << endl;
            }
        }
    }

  cout << "Total number of events = " << Total << endl;
  cout << "Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 << endl;
  cout << "Rate L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8/Total*scale_rate*scale_lumi*1000 << " Hz" << endl;

  return;

}
