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
#include <vector>
#include <TPaletteAxis.h>
#include <TLatex.h>

using namespace std;

// Read a event and returns a pointer to the list of jetpt and mjj
void ReadVBFEvent (TTree* inTree, UInt_t  i_ev, vector<Float_t>* JetPt1, vector<Float_t>* JetPt2, vector<Float_t>* Mjj, vector<Float_t>* Muon, vector<Float_t>* Tau, vector<Float_t>* MET)
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
    vector<float>   *in_l1tMuPt = 0;
    vector<float>   *in_l1tMuEta = 0;
    vector<float>   *in_l1tMuPhi = 0;
    vector<float>   *in_tauPt = 0;
    vector<float>   *in_tauEta = 0;
    vector<float>   *in_tauPhi = 0;
    vector<float>   *in_l1tauPt = 0;
    vector<float>   *in_l1tauEta = 0;
    vector<float>   *in_l1tauPhi = 0;
    Float_t         in_pfMetNoMuEt = 0;
    vector<float>   *in_l1t_MET_type = 0;
    vector<float>   *in_l1t_MET_mt = 0;


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
    inTree->SetBranchAddress("l1t_muons_pt", &in_l1tMuPt);
    inTree->SetBranchAddress("l1t_muons_eta", &in_l1tMuEta);
    inTree->SetBranchAddress("l1t_muons_phi", &in_l1tMuPhi);
    inTree->SetBranchAddress("tauPt", &in_tauPt);
    inTree->SetBranchAddress("tauEta", &in_tauEta);
    inTree->SetBranchAddress("tauPhi", &in_tauPhi);
    inTree->SetBranchAddress("l1tPt", &in_l1tauPt);
    inTree->SetBranchAddress("l1tEta", &in_l1tauEta);
    inTree->SetBranchAddress("l1tPhi", &in_l1tauPhi);
    inTree->SetBranchAddress("pfMetNoMuEt", &in_pfMetNoMuEt);
    inTree->SetBranchAddress("l1t_MET_type", &in_l1t_MET_type);
    inTree->SetBranchAddress("l1t_MET_mt", &in_l1t_MET_mt);

    inTree->GetEntry(i_ev);

    // Check if there are at least 2 jets, 1 tau and 1 mu both in the online and offline objects
    bool offline_objects = in_PxJet->size() > 1 && in_MuPt->size() > 0 && in_tauPt->size() > 0 ;
    bool L1_objects = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 ;

    if (offline_objects && L1_objects)
      {

        Float_t jet1_pt = sqrt(pow(in_PxJet->at(0),2)+pow(in_PyJet->at(0),2)+pow(in_PzJet->at(0),2));
        JetPt1->push_back(jet1_pt);
        Float_t jet2_pt = sqrt(pow(in_PxJet->at(1),2)+pow(in_PyJet->at(1),2)+pow(in_PzJet->at(1),2));
        JetPt2->push_back(jet2_pt);
        Muon->push_back(in_MuPt->at(0));
        Tau->push_back(in_tauPt->at(0));

        TLorentzVector Offline_Jet1;
        Float_t jet1_energy = sqrt(pow(in_PxJet->at(0),2)+pow(in_PyJet->at(0),2)+pow(in_PzJet->at(0),2));
        Offline_Jet1.SetPxPyPzE(in_PxJet->at(0), in_PyJet->at(0), in_PzJet->at(0), jet1_energy);
        TLorentzVector Offline_Jet2;
        Float_t jet2_energy = sqrt(pow(in_PxJet->at(1),2)+pow(in_PyJet->at(1),2)+pow(in_PzJet->at(1),2));
        Offline_Jet2.SetPxPyPzE(in_PxJet->at(1), in_PyJet->at(1), in_PzJet->at(1), jet2_energy);
        TLorentzVector Offline_DiJet;
        Offline_DiJet = Offline_Jet1 + Offline_Jet2;
        Float_t m_jj = Offline_DiJet.M();
        Mjj->push_back(m_jj);
      }

    MET->push_back(in_pfMetNoMuEt);

  }

// Plot jetPt and mjj distribution for VBF events
void PlotDistributionsRaw(TString EventSample)
  {

    // TString EventSample = "MC_MiniAOD_VBFHHTo2B2Tau_12_01_23";
    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/"+EventSample;
    TString Output = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotDistributions/"+EventSample+"_Distribution";
    system("mkdir -p "+Output);

    UInt_t events_Total = 0; 
    vector<Float_t> JetPt1;
    vector<Float_t> JetPt2;
    vector<Float_t> Mjj;
    vector<Float_t> Muon;
    vector<Float_t> Tau;
    vector<Float_t> MET;

    // loop over all events in the ntuples to check which one passes the VBF or MuTau selections
    for (int nt = 0 ; nt < 5 ; ++nt)
      {
        TString FileName_ntuples = Path_ntuples + "/Ntuple_" + to_string(nt) + ".root";
        TFile f_in(FileName_ntuples.Data(), "READ");
        TDirectoryFile* df = (TDirectoryFile*)f_in.Get("Ntuplizer_noTagAndProbe_multipleTaus");
        TTree* inTree = (TTree*)df->Get("TagAndProbe");

        cout << "\nNumber of events in the Ntuple " << nt << " = " << inTree->GetEntries() << endl;
        events_Total += inTree->GetEntries();

        for (UInt_t i_ev = 0 ; i_ev < inTree->GetEntries() ; ++i_ev)
          {
            ReadVBFEvent (inTree, i_ev, &JetPt1, &JetPt2, &Mjj, &Muon, &Tau, &MET);
          }
      }

    // gStyle->SetLabelSize(.05, "XY");
    // TLegend* Legend =  new TLegend(0.15,0.75,0.48,0.88);
    // Legend->SetBorderSize(0);
    // Legend->AddEntry(Histo_JetPt1 , "p_{T} jet1 [GeV]", "LPE");
    // Legend->Draw();

    // plot distributions JetPt1
    auto min_JetPt1 = *min_element(JetPt1.begin(), JetPt1.end());
    auto max_JetPt1 = *max_element(JetPt1.begin(), JetPt1.end());  
    TH1F* Histo_JetPt1 = new TH1F("Histo_JetPt1", "Histo_JetPt1", 50, 0, max_JetPt1);
    for (UInt_t i = 0; i < JetPt1.size(); ++i)
      {
        Histo_JetPt1->Fill(JetPt1.at(i));
      }
 
    TH1F* Histo_JetPt2 = new TH1F("Histo_JetPt2", "Histo_JetPt2", 50, 0, max_JetPt1);
    for (UInt_t i = 0; i < JetPt2.size(); ++i)
      {
        Histo_JetPt2->Fill(JetPt2.at(i));
      }

    TCanvas c1 ("c1", "c1", 900., 650.);
    c1.SetGrid(10,10);
    Histo_JetPt1->SetTitle("");
    Histo_JetPt1->SetStats(0);
    Histo_JetPt1->Draw();
    Histo_JetPt1->SetLineWidth(2);
    Histo_JetPt1->SetLineColor(kBlue);
    Histo_JetPt2->Draw("same");
    Histo_JetPt2->SetLineWidth(2);
    Histo_JetPt2->SetLineColor(kRed);
    Histo_JetPt1->GetXaxis()->SetTitle("Offline Jet p_{T} [GeV]");
    // Histo_JetPt1->GetXaxis()->SetTitleOffset(1.3);
    Histo_JetPt1->GetXaxis()->SetTitleSize(0.04);
    Histo_JetPt1->GetYaxis()->SetTitle("entries");
    Histo_JetPt1->GetYaxis()->SetTitleOffset(1.3);
    Histo_JetPt1->GetYaxis()->SetTitleSize(0.04);
    Histo_JetPt1->GetYaxis()->SetRangeUser(0,1.05*Histo_JetPt2->GetMaximum());

    TLegend* Legend =  new TLegend(0.65,0.8,0.88,0.88);
    Legend->SetBorderSize(0);
    Legend->AddEntry(Histo_JetPt1 , "Jet 1", "LPE");
    Legend->AddEntry(Histo_JetPt2 , "Jet 2", "LPE");
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

    c1.SaveAs(Output+"/Distribution_JetPt.png");
    c1.SaveAs(Output+"/Distribution_JetPt.pdf");
    c1.Close();

    // plot distributions Mjj
    auto min_Mjj = *min_element(Mjj.begin(), Mjj.end());
    auto max_Mjj = *max_element(Mjj.begin(), Mjj.end());  
    TH1F* Histo_Mjj = new TH1F("Histo_Mjj", "Histo_Mjj", 50, 0, max_Mjj);
    for (UInt_t i = 0; i < Mjj.size(); ++i)
      {
        Histo_Mjj->Fill(Mjj.at(i));
      }
    TCanvas c3 ("c3", "c3", 900., 650.);
    c3.SetGrid(10,10);
    Histo_Mjj->SetTitle("");
    Histo_Mjj->SetStats(0);
    Histo_Mjj->Draw();
    Histo_Mjj->SetLineWidth(2);
    Histo_Mjj->SetLineColor(kBlue);
    Histo_Mjj->GetXaxis()->SetTitle("M_{jj} [GeV]");
    // Histo_Mjj->GetXaxis()->SetTitleOffset(1.3);
    Histo_Mjj->GetXaxis()->SetTitleSize(0.04);
    Histo_Mjj->GetYaxis()->SetTitle("entries");
    Histo_Mjj->GetYaxis()->SetTitleOffset(1.3);
    Histo_Mjj->GetYaxis()->SetTitleSize(0.04);
    // Histo_Mjj->GetXaxis()->SetRangeUser();

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

    // plot distributions Muon
    auto min_Muon = *min_element(Muon.begin(), Muon.end());
    auto max_Muon = *max_element(Muon.begin(), Muon.end());  
    TH1F* Histo_Muon = new TH1F("Histo_Muon", "Histo_Muon", 100, 0, 500);
    for (UInt_t i = 0; i < Muon.size(); ++i)
      {
        Histo_Muon->Fill(Muon.at(i));
      }
    TCanvas c4 ("c4", "c4", 900., 650.);
    c4.SetGrid(10,10);
    Histo_Muon->SetTitle("");
    Histo_Muon->SetStats(0);
    Histo_Muon->Draw();
    Histo_Muon->SetLineWidth(2);
    Histo_Muon->SetLineColor(kBlue);
    Histo_Muon->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV]");
    // Histo_Muon->GetXaxis()->SetTitleOffset(1.3);
    Histo_Muon->GetXaxis()->SetTitleSize(0.04);
    Histo_Muon->GetYaxis()->SetTitle("entries");
    Histo_Muon->GetYaxis()->SetTitleOffset(1.3);
    Histo_Muon->GetXaxis()->SetRangeUser(min_Muon, 200);
    // Histo_Muon->GetXaxis()->SetRangeUser();

    TLatex Tex5;
    Tex5.SetTextSize(0.03);
    Tex5.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex5.Draw("same");

    TLatex Tex6;
    Tex6.SetTextSize(0.035);
    Tex6.SetTextAlign(31);
    Tex6.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex6.Draw("same");

    c4.SaveAs(Output+"/Distribution_MuonPt.png");
    c4.SaveAs(Output+"/Distribution_MuonPt.pdf");

    Histo_Muon->GetXaxis()->SetRangeUser(min_Muon, 500);
    gPad->SetLogy();
    c4.SaveAs(Output+"/Distribution_MuonPt_log.png");
    c4.SaveAs(Output+"/Distribution_MuonPt_log.pdf");
    c4.Close();

    // plot distributions Tau
    auto min_Tau = *min_element(Tau.begin(), Tau.end());
    auto max_Tau = *max_element(Tau.begin(), Tau.end());  
    TH1F* Histo_Tau = new TH1F("Histo_Tau", "Histo_Tau", 50, 0, max_Tau);
    for (UInt_t i = 0; i < Tau.size(); ++i)
      {
        Histo_Tau->Fill(Tau.at(i));
      }
    TCanvas c5 ("c5", "c5", 900., 650.);
    c5.SetGrid(10,10);
    Histo_Tau->SetTitle("");
    Histo_Tau->SetStats(0);
    Histo_Tau->Draw();
    Histo_Tau->SetLineWidth(2);
    Histo_Tau->SetLineColor(kBlue);
    Histo_Tau->GetXaxis()->SetTitle("p_{T}^{#tau} [GeV]");
    // Histo_Tau->GetXaxis()->SetTitleOffset(1.3);
    Histo_Tau->GetXaxis()->SetTitleSize(0.04);
    Histo_Tau->GetYaxis()->SetTitle("entries");
    Histo_Tau->GetYaxis()->SetTitleOffset(1.3);
    Histo_Tau->GetYaxis()->SetTitleSize(0.04);
    // Histo_Tau->GetXaxis()->SetRangeUser();

    Tex3.SetTextSize(0.03);
    Tex3.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex3.Draw("same");

    Tex4.SetTextSize(0.035);
    Tex4.SetTextAlign(31);
    Tex4.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex4.Draw("same");

    c5.SaveAs(Output+"/Distribution_TauPt.png");
    c5.SaveAs(Output+"/Distribution_TauPt.pdf");
    c5.Close();

    // plot distributions MET
    auto min_MET = *min_element(MET.begin(), MET.end());
    auto max_MET = *max_element(MET.begin(), MET.end());  
    TH1F* Histo_MET = new TH1F("Histo_MET", "Histo_MET", 50, 0, 2000);
    for (UInt_t i = 0; i < MET.size(); ++i)
      {
        Histo_MET->Fill(MET.at(i));
      }
    TCanvas c6 ("c6", "c6", 900., 650.);
    c6.SetGrid(10,10);
    Histo_MET->SetTitle("");
    Histo_MET->SetStats(0);
    Histo_MET->Draw();
    Histo_MET->SetLineWidth(2);
    Histo_MET->SetLineColor(kBlue);
    Histo_MET->GetXaxis()->SetTitle("MET [GeV]");
    // Histo_MET->GetXaxis()->SetTitleOffset(1.3);
    Histo_MET->GetXaxis()->SetTitleSize(0.04);
    Histo_MET->GetYaxis()->SetTitle("entries");
    Histo_MET->GetYaxis()->SetTitleOffset(1.3);
    Histo_MET->GetYaxis()->SetTitleSize(0.04);
    // Histo_MET->GetXaxis()->SetRangeUser();

    Tex3.SetTextSize(0.03);
    Tex3.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex3.Draw("same");

    Tex4.SetTextSize(0.035);
    Tex4.SetTextAlign(31);
    Tex4.DrawLatexNDC(0.90,0.91,"(14 TeV)");
    Tex4.Draw("same");

    c6.SaveAs(Output+"/Distribution_MET.png");
    c6.SaveAs(Output+"/Distribution_MET.pdf");
    gPad->SetLogy();

    c6.SaveAs(Output+"/Distribution_MET_log.png");
    c6.SaveAs(Output+"/Distribution_MET_log.pdf");
    c6.Close();


  }