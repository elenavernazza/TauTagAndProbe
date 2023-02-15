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
#include "../Utils/CheckL1Triggers.h"

using namespace std;

// Loops over the VBF events in MC_MiniAOD_multipleTaus_15_12_21 ntuples and computes the VBF, MuTau and common acceptance
// It returns 3D plots, one for each value of muonPt, of the acceptance of the VBF trigger as a function of the other 3 dimesnions
// It computes the VBF trigger values giving the maximum acceptance and the corresponding gain
// void PlotAcceptance(Int_t starting_x_bin)
void PlotGain()
  {

    Int_t starting_x_bin = 1;

    TString FileName_rate = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/FindGoodRates/Run362616_L1_DoubleJet_X_Y_Mass_MinZ_30_30_MuW_OpenQual/Rates_4D.root";
    cout << "\nReading rates from :\n" << FileName_rate << endl;

    TFile f (FileName_rate.Data(),"READ");
    THnF* rate_4D = (THnF*)f.Get("rates_4D");

    vector<array<Float_t, 4>> set_of_on_cuts;
    vector<array<Float_t, 4>> set_of_off_cuts;
    vector<array<Int_t, 4>> set_of_on_bins;

    // definition of online and offline selections list giving a good rate close to 1 kHz 
    FindGoodRates(rate_4D, &set_of_on_cuts, &set_of_off_cuts, &set_of_on_bins, starting_x_bin); // defined in ../Utils/CheckL1Triggers.h

    cout << "Size of the array with Online cuts = " << set_of_on_cuts.size() << endl;
    cout << "Size of the array with Offline cuts = " << set_of_off_cuts.size() << endl;

    UInt_t events_Total = 0; 
    UInt_t acceptance_MuTau = 0;
    vector<UInt_t> acceptance_VBF;
    vector<UInt_t> acceptance_MuTau_VBF;

    for (UInt_t i_cut = 0; i_cut < set_of_on_cuts.size(); ++i_cut)
      {
        acceptance_VBF.push_back(0.);
        acceptance_MuTau_VBF.push_back(0.);
      }

    TString JetIDType = "tightLepVetoJetID"; // if not specified it's _tightLepVetoJetID
    TString Method = "Way1"; // if not specified it's _Way1
    TString JetSel30 = "JetSel30"; // if not specified it's _JetSel30
    bool checkJetSel30 = false;
    if (JetSel30 == "JetSel30") {checkJetSel30 = true;}
    else if (JetSel30 == "NoJetSel30") {checkJetSel30 = false;}
    else {cout << "Chose another name! " << endl;}
    
    TString EventSample = "MC_MiniAOD_VBFHHTo2B2Tau_14_02_23";
    TString Path_ntuples = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/"+EventSample;
    TString Output = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/PlotGain_wrt_Mu18Tau24/"+EventSample+"_GainPlots";
    system("mkdir -p "+Output);

    cout << "\nReading events from :\n" << Path_ntuples << endl;
    cout << "\nSaving results in :\n" << Output << endl;

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
        // for (UInt_t i_ev = 0 ; i_ev < 20 ; ++i_ev)
          {
            // The MuTau trigger is computed separately, since it has fixed thresholds values
            bool pass_MuTau = CheckMuTau(inTree, i_ev, JetIDType, Method, checkJetSel30); // defined in ../Utils/CheckL1Triggers.h
            if (pass_MuTau)
              {
                ++ acceptance_MuTau;
              }
            // The VBF trigger is computed for all the possible good combinations of selections
            CheckVBF_vs_MuTau(inTree, i_ev, set_of_on_cuts, set_of_off_cuts, pass_MuTau, &acceptance_VBF, &acceptance_MuTau_VBF, JetIDType, Method, checkJetSel30); // defined in ../Utils/CheckL1Triggers.h

          }
      }

    Int_t bins[4] = {rate_4D->GetAxis(0)->GetNbins(), rate_4D->GetAxis(1)->GetNbins(), rate_4D->GetAxis(2)->GetNbins(), rate_4D->GetAxis(3)->GetNbins()};
    Double_t xmin[4] = {rate_4D->GetAxis(0)->GetBinLowEdge(1), rate_4D->GetAxis(1)->GetBinLowEdge(1), rate_4D->GetAxis(2)->GetBinLowEdge(1), rate_4D->GetAxis(3)->GetBinLowEdge(1)};
    Double_t xmax[4] = {rate_4D->GetAxis(0)->GetBinUpEdge(bins[0]), rate_4D->GetAxis(1)->GetBinUpEdge(bins[1]), rate_4D->GetAxis(2)->GetBinUpEdge(bins[2]), rate_4D->GetAxis(3)->GetBinUpEdge(bins[3])};

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

    TString rootname = Output+"/Gain_3D_Mu18Tau24_PtJet1_PtJet2_Mjj_PtMuon.root";
    TFile* out_file = TFile::Open(rootname, "recreate");
    out_file->cd();

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
        // PtJet1_PtJet2_Mjj.at(m)->SetMinimum(0.4);
        // PtJet1_PtJet2_Mjj.at(m)->SetMaximum(1.0);
        PtJet1_PtJet2_Mjj.at(m)->SetMinimum(0.);
        PtJet1_PtJet2_Mjj.at(m)->SetMaximum(max(1.,1.3*PtJet1_PtJet2_Mjj.at(m)->GetMaximum()));
        PtJet1_PtJet2_Mjj.at(m)->Write();

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
        Tex3.DrawLatexNDC(0.5,0.93,Form("Gain (p_{T}^{#mu} > %d [GeV]) wrt Mu18Tau24", muon_pt));
        Tex3.Draw("same");

        gPad->Update();
        TPaletteAxis* pal = (TPaletteAxis*) PtJet1_PtJet2_Mjj.at(m)->GetListOfFunctions()->FindObject("palette");
        pal->SetLabelSize(0.01);
        pal->SetX1NDC(0.916); // It doesn't work (Why??)
        // PtJet1_PtJet2_Mjj.at(m)->SaveAs(Form(Output+"/Gain_PtJet1_PtJet2_Mjj_PtMuon%i.root", muon_pt));
        c.SaveAs(Form(Output+"/Gain_PtJet1_PtJet2_Mjj_PtMuon%d.png", muon_pt));
        c.SaveAs(Form(Output+"/Gain_PtJet1_PtJet2_Mjj_PtMuon%d.pdf", muon_pt));
      }

    out_file->Close();

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

    PlotGain_2D_ptj1_mjj (30., 3., bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], set_of_on_cuts, set_of_on_bins, acceptance_VBF, acceptance_MuTau_VBF, acceptance_MuTau, Output); // defined in ../Utils/CheckL1Triggers.h
    PlotGain_2D_ptj1_mjj (30., 7., bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], set_of_on_cuts, set_of_on_bins, acceptance_VBF, acceptance_MuTau_VBF, acceptance_MuTau, Output); // defined in ../Utils/CheckL1Triggers.h

  }