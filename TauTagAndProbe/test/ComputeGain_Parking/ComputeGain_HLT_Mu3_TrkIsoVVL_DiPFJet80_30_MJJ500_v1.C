//Compile with:
//c++ `root-config --cflags --libs --evelibs` -Wno-unsequenced -O3 ComputeGain_HLT_Mu3_TrkIsoVVL_DiPFJet80_30_MJJ500_v1.cc -o simple_loop

#include <iostream>
#include <map>
#include <string>

#include <TFile.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>

using namespace std;

//Do a simple array-based loop on NanoAOD and saves the output to a TFile in out.root
//More info on NanoAOD is located here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html
void loop_plain(TTreeReader& reader) {

    reader.Restart();
    
    TStopwatch sw;
    
    const unsigned int nbytes = 0;
   
    //As usual in C++ ROOT, we have to predefine every branch and this can get annoying
    //You can find the list of all branches with description here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html
    TTreeReaderValue<unsigned int> nJet(reader, "nJet");
    TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<float> Jet_mass(reader, "Jet_mass");

    TTreeReaderValue<unsigned int> nMuon(reader, "nMuon");
    TTreeReaderArray<bool> Muon_looseId(reader, "Muon_looseId");
    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");

    TTreeReaderArray<float> MET_pt(reader, "MET_pt");

    TTreeReaderValue<unsigned int> nGenPart (reader, "nGenPart");
    TTreeReaderArray<int> GenPart_pdgId (reader, "GenPart_pdgId");
    TTreeReaderArray<int> GenPart_status (reader, "GenPart_status");
    TTreeReaderArray<float> GenPart_mass (reader, "GenPart_mass");

    unsigned int nevents = 0;
    unsigned int nSUSYevents = 0;
    TFile out("out.root", "RECREATE");

    TH1D h_Jet_pt("h_Jet_pt", "h_Jet_pt", 100, 0, 1000);
    TH1D h_Jet_phi("h_Jet_phi", "h_Jet_phi", 100, 0, 1000);
    TH1D h_Jet_eta("h_Jet_eta", "h_Jet_eta", 100, 0, 1000);

    TH1D h_Muon_pt("h_Muon_pt", "h_Muon_pt", 100, 0, 1000);
    TH1D h_Muon_phi("h_Muon_phi", "h_Muon_phi", 100, 0, 1000);
    TH1D h_Muon_eta("h_Muon_eta", "h_Muon_eta", 100, 0, 1000);

    TH1D h_MET_pt("h_MET_pt", "h_MET_pt", 100, 0, 1000);

    sw.Start();

    int passing_onlyA = 0;
    int passing_B = 0;

    while (reader.Next()) {

        nevents += 1;

        ////////////////////////// begin dM filter//////////////////////////

        //Initialize particleID & Status to 0
        int particle_id = 0;
        int particle_status = 0;
        double particle_mass = 0;

        //Choose EWKino masses
        double N1_mass = 290.0;
        double N2_mass = 300.0;
        double C1_mass = 300.0;
        double epsilon = 0.0001;

        bool SUSY_event = false;

        for (unsigned int _nGenP = 0; _nGenP < *nGenPart; _nGenP ++) {
            particle_id = GenPart_pdgId[_nGenP];
            particle_status = GenPart_status[_nGenP];
            particle_mass = GenPart_mass[_nGenP];

            //pdgID 1000023 is N2
            //pdgID 1000024 is C1
            //pdgID 1000022 is N1
            //status 62 is "outgoing subprocess particle with primordial kT included"
            //status 1 is stable
            //info here https://pythia.org/latest-manual/ParticleProperties.html
            //Particle index gets added to vector genSUSYPartIndex

            bool cond1 = (particle_id == 1000023 && particle_status == 62 && (abs(particle_mass - N2_mass) < epsilon));
            bool cond2 = (particle_id == 1000024 && particle_status == 62 && (abs(particle_mass - C1_mass) < epsilon));
            bool cond3 = (particle_id == 1000022 && particle_status == 1  && (abs(particle_mass - N1_mass) < epsilon));
            if ( cond1 | cond2 | cond3 ){
                SUSY_event = true;
                // cout << "The if condition is passed" << endl;
                break;
            }

            // cout << "Event " << nevents << " nGen " << _nGenP << endl;
            // cout << "GenPart_pdgId " << GenPart_pdgId[_nGenP] << endl;
            // cout << "GenPart_status " << GenPart_status[_nGenP] << endl;
            // cout << "GenPart_mass " << GenPart_mass[_nGenP] << endl;
        }

        if (!SUSY_event) {continue;}
        nSUSYevents += 1;

        ////////////////////////// end dM filter//////////////////////////

        h_Jet_pt.Fill(Jet_pt[0]);
        h_Jet_phi.Fill(Jet_phi[0]);
        h_Jet_eta.Fill(Jet_eta[0]);

        h_Muon_pt.Fill(Muon_pt[0]);
        h_Muon_phi.Fill(Muon_phi[0]);
        h_Muon_eta.Fill(Muon_eta[0]);

        h_MET_pt.Fill(MET_pt[0]);

        // Muons
        TLorentzVector MyGoodMuon;
        MyGoodMuon.SetPtEtaPhiM(0,0,0,0.105);
        TLorentzVector MyGoodOpenMuon;
        MyGoodOpenMuon.SetPtEtaPhiM(0,0,0,0.105);
        for (unsigned int _nMuon = 0; _nMuon < *nMuon; _nMuon ++) {
            TLorentzVector MyMuon;
            MyMuon.SetPtEtaPhiM(Muon_pt[_nMuon],Muon_eta[_nMuon],Muon_phi[_nMuon],0.105);
            if (MyMuon.Pt() > MyGoodMuon.Pt()) {
                MyGoodMuon = MyMuon;
            }
            if ((MyMuon.Pt() > MyGoodOpenMuon.Pt()) & Muon_looseId[_nMuon]) {
                MyGoodOpenMuon = MyMuon;
            }
        }

        TLorentzVector MyGoodJet1;
        MyGoodJet1.SetPtEtaPhiM(0,0,0,0);
        unsigned int _nJet1_id = 0;
        for (unsigned int _nJet1 = 0; _nJet1 < *nJet; _nJet1 ++) {
            TLorentzVector MyJet1;
            MyJet1.SetPtEtaPhiM(Jet_pt[_nJet1],Jet_eta[_nJet1],Jet_phi[_nJet1],0.);
            if (MyJet1.Pt() > MyGoodJet1.Pt()) {
                MyGoodJet1 = MyJet1;
                _nJet1_id = _nJet1;
            }
        }
        TLorentzVector MyGoodJet2;
        MyGoodJet2.SetPtEtaPhiM(0,0,0,0);
        for (unsigned int _nJet2 = 0; _nJet2 < *nJet; _nJet2 ++) {
            if (_nJet2 == _nJet1_id) continue;
            TLorentzVector MyJet2;
            MyJet2.SetPtEtaPhiM(Jet_pt[_nJet2],Jet_eta[_nJet2],Jet_phi[_nJet2],0.);
            if (MyJet2.Pt() > MyGoodJet2.Pt()) {
                MyGoodJet2 = MyJet2;
            } 
        }

        TLorentzVector MyGoodDiJet40;
        MyGoodDiJet40.SetPtEtaPhiM(0,0,0,0);
        for (unsigned int _nJet1 = 0; _nJet1 < *nJet; _nJet1 ++) {
            if (Jet_pt[_nJet1] < 40) continue;
            TLorentzVector MyJet1_40;
            MyJet1_40.SetPtEtaPhiM(Jet_pt[_nJet1],Jet_eta[_nJet1],Jet_phi[_nJet1],0.);
            for (unsigned int _nJet2 = _nJet1+1; _nJet2 < *nJet; _nJet2 ++) {
                if (Jet_pt[_nJet2] < 40) continue;
                TLorentzVector MyJet2_40;
                MyJet2_40.SetPtEtaPhiM(Jet_pt[_nJet2],Jet_eta[_nJet2],Jet_phi[_nJet2],0.);
                TLorentzVector MyDiJet40;
                MyDiJet40 = MyJet1_40 + MyJet2_40;
                if (MyDiJet40.M() > MyGoodDiJet40.M()) {
                    MyGoodDiJet40 = MyDiJet40;
                }
            }
        }
        TLorentzVector MyGoodDiJet45;
        MyGoodDiJet45.SetPtEtaPhiM(0,0,0,0);
        for (unsigned int _nJet1 = 0; _nJet1 < *nJet; _nJet1 ++) {
            if (Jet_pt[_nJet1] < 45) continue;
            TLorentzVector MyJet1_45;
            MyJet1_45.SetPtEtaPhiM(Jet_pt[_nJet1],Jet_eta[_nJet1],Jet_phi[_nJet1],0.);
            for (unsigned int _nJet2 = _nJet1+1; _nJet2 < *nJet; _nJet2 ++) {
                if (Jet_pt[_nJet2] < 45) continue;
                TLorentzVector MyJet2_45;
                MyJet2_45.SetPtEtaPhiM(Jet_pt[_nJet2],Jet_eta[_nJet2],Jet_phi[_nJet2],0.);
                TLorentzVector MyDiJet45;
                MyDiJet45 = MyJet1_45 + MyJet2_45;
                if (MyDiJet45.M() > MyGoodDiJet45.M()) {
                    MyGoodDiJet45 = MyDiJet45;
                }
            }
        }


        // New parking strategy
        // HLT_VBFparking_Mu3_TrkIsoVVL_DiPFJet80_30_Mjj500_v
        bool A_muon_open_3 = MyGoodOpenMuon.Pt() > 3;
        bool A_jet1_100 = MyGoodJet1.Pt() > 100;
        bool A_jet2_40 = MyGoodJet2.Pt() > 40;
        bool A_mjj_600 = MyGoodDiJet40.M() > 600;
        bool A_MET_100 = MET_pt[0] > 100;
        bool A = A_muon_open_3 & A_jet1_100 & A_jet2_40 & A_mjj_600 & A_MET_100;

        // Old MET and VBF strategies
        // 1) MET_MHT
        bool B1_met_230 = MET_pt[0] > 230;
        bool B1_muon_3 = MyGoodMuon.Pt() > 3; 
        bool B1 = B1_met_230 & B1_muon_3;
        // 2) HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60
        bool B2_muon_8 = MyGoodMuon.Pt() > 8;
        bool B2_jet1_135 = MyGoodJet1.Pt() > 135;
        bool B2_jet2_45 = MyGoodJet2.Pt() > 45;
        bool B2_mjj_850 = MyGoodDiJet45.M() > 850;
        bool B2_htt_400 = true; //[FIXME]
        bool B2_METNoMu_100 = MET_pt[0] > 100; //[FIXME]
        bool B2 = B2_muon_8 & B2_jet1_135 & B2_jet2_45 & B2_mjj_850 & B2_htt_400 & B2_METNoMu_100;
        // 3) HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v
        bool B3_muon_4 = MyGoodMuon.Pt() > 4;
        bool B3_jet1_135 = MyGoodJet1.Pt() > 135;
        bool B3_jet2_45 = MyGoodJet2.Pt() > 45;
        bool B3_mjj_850 = MyGoodDiJet45.M() > 850;
        bool B3_htt_400 = true; //[FIXME]
        bool B3_METNoMu_100 = MET_pt[0] > 100; //[FIXME]
        bool B3 = B3_muon_4 & B3_jet1_135 & B3_jet2_45 & B3_mjj_850 & B3_htt_400 & B3_METNoMu_100;
        bool B = B1 | B2 | B3;

        if (B) {
            passing_B += 1;
        }
        if (A & !B) {
            passing_onlyA += 1;
        }
        
    }

    sw.Stop();
    const auto cpu_time = sw.CpuTime();
    const auto real_time = sw.RealTime();

    cout << "\nNumber of events = " << nevents << "\nNumber of SUSY events = " << nSUSYevents << "\n" << endl;

    const auto speed = (double)nevents / real_time / 1000.0;
    cout << "h_Jet_pt " << h_Jet_pt.GetEntries() << " " << h_Jet_pt.GetMean() << endl;
    cout << "h_Muon_pt " << h_Muon_pt.GetEntries() << " " << h_Muon_pt.GetMean() << endl;
    cout << "h_MET_pt " << h_MET_pt.GetEntries() << " " << h_MET_pt.GetMean() << endl;

    cout << "\npassing_onlyA = " << passing_onlyA << endl;
    cout << "passing_B = " << passing_B << endl;
    cout << "Gain = " << float(float(passing_onlyA)/float(passing_B)) << endl;

    cout << "\nloop_plain " << "cpu_time=" << cpu_time << ",real_time=" << real_time << ",speed=" << speed << " kHz" << endl;

    out.Write();
    out.Close();
}

void Loop ()
{
    TString FileName = "387ED486-1915-1740-A901-EA0C62954DDD.root";

    gROOT->SetBatch(true);
    TFile tf(FileName);
    TTreeReader reader("Events", &tf);
    loop_plain(reader);
}