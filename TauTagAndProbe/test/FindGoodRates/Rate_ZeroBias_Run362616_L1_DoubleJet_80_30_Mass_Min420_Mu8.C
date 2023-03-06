#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
// #include "Build_Isolation_WPs.C"

// Header file for the classes stored in the TTree if any.
#include "interface/L1AnalysisL1UpgradeDataFormat.h"
// #include "interface/L1AnalysisL1CaloClusterDataFormat.h"
#include "interface/L1AnalysisEventDataFormat.h"

using namespace std;

void ReadFiredEvents (vector<ULong_t>* v_events, vector<double>* v_ptjet1, vector<double>* v_ptjet2, vector<double>* v_mjj, vector<double>* v_ptmu)
{
    TString FiredL1Events365File = "/afs/cern.ch/user/e/evernazz/L1Studies/RateComputation/CMSSW_13_0_0_pre2/src/L1Trigger/L1TNtuples/FiredL1Events365.txt";
    ifstream indata;
    indata.open(FiredL1Events365File);
    int event; double ptjet1; double ptjet2; double mjj; double ptmu;

    cout << "Start reading file " << FiredL1Events365File << endl;

    while (!indata.eof()) 
    {
        indata >> event >> ptjet1 >> ptjet2 >> mjj >> ptmu;
        v_events->push_back(event);
        v_ptjet1->push_back(ptjet1);
        v_ptjet2->push_back(ptjet2);
        v_mjj->push_back(mjj);
        v_ptmu->push_back(ptmu);
    }

    cout << "End reading" << endl;

}

void Rate ()
{

    vector<ULong_t> v_events;
    vector<double> v_ptjet1;
    vector<double> v_ptjet2;
    vector<double> v_mjj;
    vector<double> v_ptmu;

    ReadFiredEvents(&v_events, &v_ptjet1, &v_ptjet2, &v_mjj, &v_ptmu);
    bool match = false;

    Float_t Max_Events = -1;
    Float_t Total = 0;
    Float_t Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = 0.;

    cout << "Begin loop" << endl;

    TString path = "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/elfontan/condor/2022EphZB_run362616_126X";

    float nb = 2450.;
    float scale_rate = 0.001*(nb*11245.6);
    float thisLumiRun = 2.050E34;
    float scaleToLumi = 2.000E34;
    float scale_lumi = scaleToLumi/thisLumiRun;

    TString NameL1EventTree = "l1EventTree/L1EventTree";
    TString NameL1UpgradeTree = "l1UpgradeTree/L1UpgradeTree";

    TChain dataL1Event(NameL1EventTree.Data());
    TChain dataL1Upgrade(NameL1UpgradeTree.Data());

    TString FileName_in = path + "/*.root";
    dataL1Event.Add(FileName_in.Data());
    dataL1Upgrade.Add(FileName_in.Data());

    dataL1Event.SetMakeClass(1);
    dataL1Upgrade.SetMakeClass(1);

    ULong_t         event;
    UInt_t          in_l1tnJets;
    vector<float>   in_l1tPtJet;
    vector<float>   in_l1tEtaJet;
    vector<float>   in_l1tPhiJet;
    vector<float>   in_l1tMuPt;
    vector<float>   in_l1tMuEta;
    vector<float>   in_l1tMuPhi;
    vector<unsigned short int>   in_l1tMuQual;

    TBranch *b_L1Event_event;
    TBranch *b_L1Upgrade_l1tnJets;
    TBranch *b_L1Upgrade_l1tPtJet;
    TBranch *b_L1Upgrade_l1tPhiJet;
    TBranch *b_L1Upgrade_l1tEtaJet;
    TBranch *b_L1Upgrade_l1tMuPt;
    TBranch *b_L1Upgrade_l1tMuEta;
    TBranch *b_L1Upgrade_l1tMuPhi;
    TBranch *b_L1Upgrade_l1tMuQual;

    dataL1Event.SetBranchAddress("event",      &event,         &b_L1Event_event);
    dataL1Upgrade.SetBranchAddress("nJets",    &in_l1tnJets,   &b_L1Upgrade_l1tnJets);
    dataL1Upgrade.SetBranchAddress("jetEt",    &in_l1tPtJet,   &b_L1Upgrade_l1tPtJet);
    dataL1Upgrade.SetBranchAddress("jetEta",   &in_l1tEtaJet,  &b_L1Upgrade_l1tEtaJet);
    dataL1Upgrade.SetBranchAddress("jetPhi",   &in_l1tPhiJet,  &b_L1Upgrade_l1tPhiJet);
    dataL1Upgrade.SetBranchAddress("muonEt",   &in_l1tMuPt,    &b_L1Upgrade_l1tMuPt);
    dataL1Upgrade.SetBranchAddress("muonEta",  &in_l1tMuEta,   &b_L1Upgrade_l1tMuEta);
    dataL1Upgrade.SetBranchAddress("muonPhi",  &in_l1tMuPhi,   &b_L1Upgrade_l1tMuPhi);
    dataL1Upgrade.SetBranchAddress("muonQual", &in_l1tMuQual,  &b_L1Upgrade_l1tMuQual);

    int nEntries = dataL1Upgrade.GetEntries();
    cout << "TChain has " << nEntries << " events" << endl;
    if (Max_Events == -1)
    {
        Max_Events = nEntries;
    }
    cout << "Analysing " << Max_Events << " events\n" << endl;

    for(int i = 0 ; i < Max_Events ; ++i)
    {
        dataL1Event.GetEntry(i);
        dataL1Upgrade.GetEntry(i);
        ++ Total;

        match = false;
        for (UInt_t i_fired = 0; i_fired < v_events.size(); ++i_fired)
        {
            if (event == v_events.at(i_fired))
            {
                cout << v_events.at(i_fired) << " Official Recipe : ";
                cout << v_ptjet1.at(i_fired) << " " << v_ptjet2.at(i_fired) << " ";
                cout << v_mjj.at(i_fired) << " " << v_ptmu.at(i_fired) << " (" << i_fired << ")" << endl;
                match = true;
                break;
            }
        }

        bool L1_objects_VBF = in_l1tPtJet.size() > 1 && in_l1tMuPt.size() > 0 ;

        TLorentzVector myGoodOnlineMuon;
        TLorentzVector myGoodOnlineJet1;
        TLorentzVector myGoodOnlineJet2;
        TLorentzVector myGoodOnlineDiJet;
        float myGoodOnlineMjj = -1.;

        if (L1_objects_VBF)
        {

            // take the first (most enegrgetic) muon passing the open quality requirement
	    myGoodOnlineMuon.SetPtEtaPhiM(0,0,0,0);
	    for (UInt_t iMuon = 0; iMuon < in_l1tMuPt.size(); ++iMuon)
            {
                if (in_l1tMuQual.at(iMuon) < 12) continue; // single quality: 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
		TLorentzVector myOnlineMuon;
                myOnlineMuon.SetPtEtaPhiM(in_l1tMuPt.at(iMuon),in_l1tMuEta.at(iMuon),in_l1tMuPhi.at(iMuon),0.105);
		if (myOnlineMuon.Pt() > myGoodOnlineMuon.Pt())
		{
		    myGoodOnlineMuon = myOnlineMuon;
		}
            }

            // loop among the pairs of jets (80,30) giving the highest mjj
            for (UInt_t iL1Jet1 = 0 ; iL1Jet1 < in_l1tPtJet.size() ; ++iL1Jet1)
            {
                TLorentzVector myOnlineJet1;
                myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet.at(iL1Jet1),in_l1tEtaJet.at(iL1Jet1),in_l1tPhiJet.at(iL1Jet1),0.);
                if (myOnlineJet1.Pt() < 80) continue;

                for (UInt_t jL1Jet2 = iL1Jet1 + 1 ; jL1Jet2 < in_l1tPtJet.size() ; ++jL1Jet2)
                {
                    TLorentzVector myOnlineJet2;
                    myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet.at(jL1Jet2),in_l1tEtaJet.at(jL1Jet2),in_l1tPhiJet.at(jL1Jet2),0.);
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
                cout << event << " My Recipe       : ";
                cout << myGoodOnlineJet1.Pt() << " " << myGoodOnlineJet2.Pt() << " ";
                cout << myGoodOnlineDiJet.M() << " " << myGoodOnlineMuon.Pt() << " ";
                cout << in_l1tMuQual.at(0) << endl;
            }
            else
            {
                if (match)
                {
                cout << event << " My Recipe       : ";
                cout << myGoodOnlineJet1.Pt() << " " << myGoodOnlineJet2.Pt() << " ";
                cout << myGoodOnlineDiJet.M() << " " << myGoodOnlineMuon.Pt() << " ";
                cout << in_l1tMuQual.at(0) << " Not Passed" << endl;                    
                }
            }
        }

    }

  cout << "Total number of events = " << Total << endl;
  cout << "Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8 << endl;
  cout << "Rate L1_DoubleJet_80_30_Mass_Min420_Mu8 = " << Events_L1_DoubleJet_80_30_Mass_Min420_Mu8/Total*scale_rate*scale_lumi*1000;
  cout << " +/- " << sqrt(Events_L1_DoubleJet_80_30_Mass_Min420_Mu8)/Total*scale_rate*scale_lumi*1000 << " Hz" << endl;

}
