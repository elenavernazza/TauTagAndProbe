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
#include <TString.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <cstdlib>

// I want to add more statistics, but we have to be careful and have a homogeneus luminosity
// for the ZeroBias sample 321755 we can look at this page
// https://cmsoms.cern.ch/cms/runs/lumisection?cms_run=321755&cms_run_sequence=GLOBAL-RUN
// if(in_lumi<63 || in_lumi>250) thisLumiRun = (1.686 + 1.522)/2 = 1.604 → 1.600E34;
// if(in_lumi<250 || in_lumi>400) thisLumiRun = (1.522 + 1.398)/2 = 1.460 → 1.600E34;
// if(in_lumi<400 || in_lumi>600) thisLumiRun = (1.398 + 1.252)/2 = 1.375 → 1.600E34;
// if(in_lumi<600 || in_lumi>728) thisLumiRun = (1.252 + 1.165)/2 = 1.208 → 1.600E34;

using namespace std;

void PlotRate_2D_ptj1_mjj(Int_t fixed_bin_y, Int_t fixed_bin_w, THnF* rates_4D, TString output, Int_t bin0, Int_t xmin0, Int_t xmax0, Int_t bin2, Int_t xmin2, Int_t xmax2, Int_t xmin3) 
  {
    // Compute 2D_rate with a fixed muonPt value and jetPt2 value 
    TH2F* rate_2D_ptj1_mjj = new TH2F("rate_2D_ptj1_mjj","rate_2D_ptj1_mjj", bin0, xmin0, xmax0, bin2, xmin2, xmax2);
    Int_t fixed_cut_y = rates_4D->GetAxis(1)->GetBinLowEdge(fixed_bin_y);
    Int_t fixed_cut_w = xmin3 + fixed_bin_w - 1;

    Float_t rate_value = 0;
    Int_t vec_bin[4] = {0, 0, 0, 0};
    for (Int_t I_bin_x = 1 ; I_bin_x <= bin0 ; ++I_bin_x)
      {
        for (Int_t I_bin_z = 1 ; I_bin_z <= bin2 ; ++I_bin_z)
          {
            Int_t vec_bin[4] = {I_bin_x, fixed_bin_y, I_bin_z, fixed_bin_w};
            rate_value = rates_4D->GetBinContent(vec_bin);
            rate_2D_ptj1_mjj->SetBinContent(I_bin_x, I_bin_z, rate_value);
          }
      }

    TCanvas* c_2D_1 = new TCanvas("c_2D_1","c_2D_1",700.,550.);
    c_2D_1->cd();
    c_2D_1->SetRightMargin(0.16); // important for labels going outside the canvas!
    rate_2D_ptj1_mjj->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
    rate_2D_ptj1_mjj->GetXaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->GetYaxis()->SetTitle("M_{jj} [GeV]");
    rate_2D_ptj1_mjj->GetYaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->GetZaxis()->SetTitle("rate [kHz]");
    rate_2D_ptj1_mjj->GetZaxis()->SetTitleOffset(1.3);
    rate_2D_ptj1_mjj->SetTitle("");
    rate_2D_ptj1_mjj->SetStats(0);
    rate_2D_ptj1_mjj->Draw("COLZ");
    rate_2D_ptj1_mjj->SetMinimum(0);
    rate_2D_ptj1_mjj->SetMaximum(3);

    TLatex Tex11;
    Tex11.SetTextSize(0.03);
    Tex11.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
    Tex11.Draw("same");

    TLatex Tex22;
    Tex22.SetTextSize(0.035);
    Tex22.SetTextAlign(31);
    Tex22.DrawLatexNDC(0.88,0.91,"(14 TeV)");
    Tex22.Draw("same");

    TLatex Tex33;
    Tex33.SetTextSize(0.035);
    Tex33.SetTextAlign(12);
    TString printout = "Zero Bias L1 rate p_{T}^{jet2} > " + to_string(fixed_cut_y) + " GeV && p_{T}^{#mu} > " + to_string(fixed_cut_w) + " GeV";
    Tex33.DrawLatexNDC(0.25,0.96, printout);
    Tex33.Draw("same");

    c_2D_1->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_" + to_string(fixed_cut_y) + "_mjj_Y_ptmu_" + to_string(fixed_cut_w) + ".png");
    c_2D_1->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_" + to_string(fixed_cut_y) + "_mjj_Y_ptmu_" + to_string(fixed_cut_w) + ".pdf");

    TH2F *clone_rate_2D_ptj1_mjj = (TH2F*) rate_2D_ptj1_mjj->Clone();
    clone_rate_2D_ptj1_mjj->SetName("clone_rate_2D_ptj1_mjj");
    // draw contour line
    double contours[2];
    contours[0] = 0.95;
    contours[1] = 1.05;
    clone_rate_2D_ptj1_mjj->SetContour(2,contours);
    clone_rate_2D_ptj1_mjj->Draw("cont3 list same");
    clone_rate_2D_ptj1_mjj->SetLineColor(kRed);
    clone_rate_2D_ptj1_mjj->SetLineStyle(kRed);

    c_2D_1->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_" + to_string(fixed_cut_y) + "_mjj_Y_ptmu_" + to_string(fixed_cut_w) + "_line.png");
    c_2D_1->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_" + to_string(fixed_cut_y) + "_mjj_Y_ptmu_" + to_string(fixed_cut_w) + "_line.pdf");
    c_2D_1->Close();
  }


void PureRate()
{

  Int_t bins[4] = {14, 14, 30, 12};
  Double_t xmin[4] = {30., 30., 200., 3.};
  Double_t xmax[4] = {100., 100., 800., 15.};

  Int_t Denominator = 0;
  UInt_t noObjects_events = 0.;
  UInt_t GoodEvents = 0.;
  vector <UInt_t> GoodEvents_Mu (bins[3], 0);
  UInt_t noLumi_events = 0.;

  cout << "Begin loop" << endl;

  TString path = "/grid_mnt/data__data.polcms/cms/vernazza/Ntuples/ZeroBias_Run321755/";
  TString output = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_2_1/src/TauTagAndProbe/TauTagAndProbe/test/FindGoodRates/PureRatesStudies";
  system("mkdir -p "+output);

  vector <TH3F*> ptjet1_ptjet2_jetmass_ptmu;
  TH2F* PtJet1_PtJet2 = new TH2F("PtJet1_PtJet2", "PtJet1_PtJet2", bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1]);

  for (Int_t m =0; m < bins[3]; m++)
    {
      ptjet1_ptjet2_jetmass_ptmu.push_back(nullptr);
    }

  for (Int_t m =0; m < bins[3]; m++)
    {
      int num = xmin[3] + m;
      TString histname = Form("ptjet1_ptjet2_jetmass_ptmu_%d", num);
      ptjet1_ptjet2_jetmass_ptmu.at(m) = new TH3F(histname, histname, bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1], bins[2], xmin[2], xmax[2]);
    }

  int number_of_not_isolate_muons = 0;
  float nb = 2544.;
  float scale_rate = 0.001*(nb*11245.6);
  float thisLumiRun = 0; // this value will change according to the LS in order to correctly weight the events
  float scaleToLumi = 2.00E34;
  float scale_lumi = 0;

  for(int nt = 0 ; nt < 200 ; ++nt)
  // for(int nt = 0 ; nt < 1 ; ++nt)
    {
      string ntuple = "Ntuple_" + to_string(nt) + ".root";
      TString FileName_in = path + ntuple ;

      // cout << FileName_in << endl;
      // TString FileName_in = path + "Ntuple_0.root";

      TFile f_in(FileName_in.Data(),"READ");
      TDirectoryFile* df = (TDirectoryFile*)f_in.Get("ZeroBias");
      TTree* inTree = (TTree*)df->Get("ZeroBias");

      ULong64_t       in_EventNumber =  0;
      Int_t           in_RunNumber =  0;
      Int_t           in_lumi =  0;
      vector<float>   *in_l1tPtJet = 0;
      vector<float>   *in_l1tEtaJet = 0;
      vector<float>   *in_l1tPhiJet = 0;
      vector<float>   *in_l1tMuPt = 0;
      vector<float>   *in_l1tMuEta = 0;
      vector<float>   *in_l1tMuPhi = 0;
      vector<float>   *in_l1tauPt = 0;
      vector<float>   *in_l1tauEta = 0;
      vector<float>   *in_l1tauPhi = 0;

      inTree->SetBranchAddress("EventNumber", &in_EventNumber);
      inTree->SetBranchAddress("RunNumber", &in_RunNumber);
      inTree->SetBranchAddress("lumi", &in_lumi);
      inTree->SetBranchAddress("l1tPtJet", &in_l1tPtJet);
      inTree->SetBranchAddress("l1tEtaJet", &in_l1tEtaJet);
      inTree->SetBranchAddress("l1tPhiJet", &in_l1tPhiJet);
      inTree->SetBranchAddress("l1tMuPt", &in_l1tMuPt);
      inTree->SetBranchAddress("l1tMuEta", &in_l1tMuEta);
      inTree->SetBranchAddress("l1tMuPhi", &in_l1tMuPhi);
      inTree->SetBranchAddress("l1tPt", &in_l1tauPt);
      inTree->SetBranchAddress("l1tEta", &in_l1tauEta);
      inTree->SetBranchAddress("l1tPhi", &in_l1tauPhi);

      if (nt%10 == 0)
        {
          cout << "\nNumber of events in the Ntuple " << nt << " = " << inTree->GetEntries() << endl;
        }

      for(UInt_t i = 0 ; i < inTree->GetEntries() ; ++i)
        {
          inTree->GetEntry(i);
          // if(i%10000==0) cout << "Entry #" << i << endl; 
          if(in_lumi<63 || in_lumi>=729)
            {
              ++ noLumi_events;
              continue;
            }

          else if(in_lumi >= 63  || in_lumi < 250) thisLumiRun = 1.600E34;
          else if(in_lumi >= 250 || in_lumi < 400) thisLumiRun = 1.460E34;
          else if(in_lumi >= 400 || in_lumi < 600) thisLumiRun = 1.375E34;
          else if(in_lumi >= 600 || in_lumi < 729) thisLumiRun = 1.208E34;

          scale_lumi = scaleToLumi/thisLumiRun;

          Float_t weight = 1.;

          ++Denominator;

          bool L1_objects_VBF = in_l1tPtJet->size() > 1 && in_l1tMuPt->size() > 0 ;
          bool L1_Objects_MuTau = in_l1tauPt->size() > 0 && in_l1tMuPt->size() > 0 ;
          bool L1_Objects_SingleMu = in_l1tMuPt->size() > 0 ;

          TLorentzVector myGoodOnlineMuon;
          TLorentzVector myGoodOnlineTau;
          TLorentzVector myGoodOnlineJet1;
          TLorentzVector myGoodOnlineJet2;
          TLorentzVector myGoodOnlineDiJet;

          bool IsolatedMuon = false;

          float myGoodOnlineMjj = -1.;

          if (L1_objects_VBF)
            {

              ++ GoodEvents;

              for (UInt_t iL1Muon = 0 ; iL1Muon < in_l1tMuPt->size() ; ++iL1Muon)
                {
                  if (IsolatedMuon) break;
                  TLorentzVector myOnlineMuon;
                  myOnlineMuon.SetPtEtaPhiM(in_l1tMuPt->at(iL1Muon), in_l1tMuEta->at(iL1Muon), in_l1tMuPhi->at(iL1Muon), 0.105);
                  bool myIsolationCondition = true;
                  for (UInt_t iL1Jet = 0 ; iL1Jet < in_l1tPtJet->size() ; ++iL1Jet)
                    {
                      if (in_l1tPtJet->at(iL1Jet) > 5)
                        {
                          TLorentzVector myOnlineJet;
                          myOnlineJet.SetPtEtaPhiM(in_l1tPtJet->at(iL1Jet),in_l1tEtaJet->at(iL1Jet),in_l1tPhiJet->at(iL1Jet),0.);
                          if (myOnlineMuon.DeltaR(myOnlineJet) < 0.3)
                            {
                              // cout << "Muon matched to " << iL1Jet << " DeltaR = " << myOnlineMuon.DeltaR(myOnlineJet) << endl;
                              myIsolationCondition = false;
                            }
                        }
                    }
                  if (myIsolationCondition)
                    {
                      myGoodOnlineMuon = myOnlineMuon;
                      IsolatedMuon = true;
                      break;
                    }
                }

              if (L1_Objects_MuTau) myGoodOnlineTau.SetPtEtaPhiM(in_l1tauPt->at(0), in_l1tauEta->at(0), in_l1tauPhi->at(0), 1.776);

              for (UInt_t iL1Jet1 = 0 ; iL1Jet1 < in_l1tPtJet->size() ; ++iL1Jet1)
                {
                  TLorentzVector myOnlineJet1;
                  myOnlineJet1.SetPtEtaPhiM(in_l1tPtJet->at(iL1Jet1),in_l1tEtaJet->at(iL1Jet1),in_l1tPhiJet->at(iL1Jet1),0.);

                  for (UInt_t jL1Jet2 = iL1Jet1 + 1 ; jL1Jet2 < in_l1tPtJet->size() ; ++jL1Jet2)
                    {
                      TLorentzVector myOnlineJet2;
                      myOnlineJet2.SetPtEtaPhiM(in_l1tPtJet->at(jL1Jet2),in_l1tEtaJet->at(jL1Jet2),in_l1tPhiJet->at(jL1Jet2),0.);
                      TLorentzVector myOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                      float myOnlineMjj = myOnlineDiJet.M();

                      // bool check = in_l1tPtJet->at(iL1Jet) > 20. && in_l1tPtJet->at(jL1Jet) > 20. && jL1Jet != iL1Jet ;
                      if (myOnlineMjj > myGoodOnlineMjj)
                        {
                          myGoodOnlineMjj = myOnlineMjj;
                          myGoodOnlineJet1 = myOnlineJet1;
                          myGoodOnlineJet2 = myOnlineJet2;
                          myGoodOnlineDiJet = myOnlineJet1 + myOnlineJet2;
                        }
                    }  
                }

              PtJet1_PtJet2->Fill(myGoodOnlineJet1.Pt(),myGoodOnlineJet2.Pt(), scale_lumi);

              bool check_MuTau = L1_Objects_MuTau && (myGoodOnlineMuon.Pt() > 18) && (myGoodOnlineTau.Pt() > 24);
              bool check_SingleMu = L1_Objects_SingleMu && (myGoodOnlineMuon.Pt() > 22);

              if (!IsolatedMuon) ++number_of_not_isolate_muons;

              if (!check_MuTau && !check_SingleMu && IsolatedMuon)
                {

                  for (Int_t n_cut = 0; n_cut < bins[3]; n_cut++)
                    {
                      int muon_cut_min = n_cut + xmin[3];
                      int muon_cut_max = n_cut + 1 + xmin[3];
                      if (n_cut == 0)
                        {
                          if (myGoodOnlineMuon.Pt() < muon_cut_max)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                      else if (n_cut == bins[3]-1)
                        {
                          if (myGoodOnlineMuon.Pt() >= muon_cut_min)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                      else
                        {
                          if (myGoodOnlineMuon.Pt() >= muon_cut_min && myGoodOnlineMuon.Pt() < muon_cut_max)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(myGoodOnlineJet1.Pt(), myGoodOnlineJet2.Pt(), myGoodOnlineDiJet.M(), scale_lumi);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                    }
                }
            }

          else
            {
              ++ noObjects_events;
            }
        }
    }

  cout << "\nNumber of not isolated muons = " << number_of_not_isolate_muons << endl;

  cout << "Denominator = " << Denominator << endl;
  cout << "Number of events without jets or muons = " << noObjects_events << endl;
  cout << "Number of events out of lumi = " << noLumi_events << endl;
  cout << "Number of good events = " << GoodEvents << endl;

  // Compute 4D histogram, where at each combination of jetPt1, jetPt2, Mjj and muonPt corresponds a rate

  THnF* rates_4D = new THnF("rates_4D","rates_4D", 4, bins, xmin, xmax);
  float integral_xyzw = -1.;
  int combinations = 0;

  for (Int_t i_bin_x = 1 ; i_bin_x <= bins[0] ; ++i_bin_x)
    {
      for (Int_t i_bin_y = 1 ; i_bin_y <= bins[1] ; ++i_bin_y)
        {
          for (Int_t i_bin_z = 1 ; i_bin_z <= bins[2] ; ++i_bin_z)
            {
              vector <float> int_xyz (bins[3], 0);
              for (Int_t w = 0; w < bins[3] ; w++)
                {
                  int_xyz.at(w) = ptjet1_ptjet2_jetmass_ptmu.at(w)->Integral(i_bin_x,bins[0]+1,i_bin_y,bins[1]+1,i_bin_z,bins[2]+1)/Denominator*scale_rate;
                }

              for (Int_t i_bin_w = 1 ; i_bin_w <= bins[3] ; ++i_bin_w)
                {
                  integral_xyzw = 0;
                  for (Int_t i_vec_w = 0; i_vec_w < bins[3]; ++i_vec_w)
                    {
                      if (i_vec_w >= i_bin_w - 1)
                        {
                          integral_xyzw += int_xyz.at(i_vec_w);
                        }
                    }

                  if (integral_xyzw < 1.05 && integral_xyzw > 0.95 && i_bin_x >= i_bin_y)
                    {
                      ++ combinations;
                    }

                  Int_t v_bin[4] = {i_bin_x, i_bin_y, i_bin_z, i_bin_w};
                  rates_4D->SetBinContent(v_bin, integral_xyzw);

                }
            }
        }
    }

  // --------------- Plotting ---------------

  PlotRate_2D_ptj1_mjj(1, 1, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(1, 2, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(1, 3, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(1, 4, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(1, 5, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(1, 6, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);
  PlotRate_2D_ptj1_mjj(3, 6, rates_4D, output, bins[0], xmin[0], xmax[0], bins[2], xmin[2], xmax[2], xmin[3]);

  // Compute 2D_rate with a fixed muonPt value and jetPt2 value 
  TH2F* rate_2D = new TH2F("rate_2D","rate_2D", bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1]);

  Int_t fixed_bin_z = 10;
  Int_t fixed_cut_z = rates_4D->GetAxis(2)->GetBinLowEdge(fixed_bin_z);
  Int_t fixed_bin_w = 6;
  Int_t fixed_cut_w = xmin[3] + fixed_bin_w - 1;

  for (Int_t I_bin_x = 1 ; I_bin_x <= bins[0] ; ++I_bin_x)
    {
      for (Int_t I_bin_y = 1 ; I_bin_y <= bins[1] ; ++I_bin_y)
        {
          Int_t vec_bin[4] = {I_bin_x, I_bin_y, fixed_bin_z, fixed_bin_w};
          Float_t rate_value = rates_4D->GetBinContent(vec_bin);
          rate_2D->SetBinContent(I_bin_x, I_bin_y, rate_value);
        }
    }

  TCanvas* c_2D = new TCanvas("c_2D","c_2D",700.,550.);
  c_2D->cd();
  c_2D->SetRightMargin(0.16); // important for labels going outside the canvas!
  rate_2D->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
  rate_2D->GetXaxis()->SetTitleOffset(1.3);
  rate_2D->GetYaxis()->SetTitle("p_{T}^{j2} [GeV]");
  rate_2D->GetYaxis()->SetTitleOffset(1.3);
  rate_2D->GetZaxis()->SetTitle("rate [kHz]");
  rate_2D->GetZaxis()->SetTitleOffset(1.3);
  rate_2D->SetTitle("");
  rate_2D->SetMaximum(2);
  rate_2D->SetStats(0);
  rate_2D->Draw("COLZ");

  TLatex Tex1;
  Tex1.SetTextSize(0.03);
  Tex1.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
  Tex1.Draw("same");

  TLatex Tex2;
  Tex2.SetTextSize(0.035);
  Tex2.SetTextAlign(31);
  Tex2.DrawLatexNDC(0.88,0.91,"(14 TeV)");
  Tex2.Draw("same");

  TLatex Tex3;
  Tex3.SetTextSize(0.035);
  Tex3.SetTextAlign(12);
  TString printout = "Zero Bias L1 rate M_{jj} > " + to_string(fixed_cut_z) + " GeV && p_{T}^{#mu} > " + to_string(fixed_cut_w) + " GeV";
  Tex3.DrawLatexNDC(0.25,0.96,printout);
  Tex3.Draw("same");

  c_2D->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_Y_mjj_" + to_string(fixed_cut_z) + "_ptmu_" + to_string(fixed_cut_w) + ".png");
  c_2D->SaveAs(output+"/pure_rate_2D_ptj1_X_ptj2_Y_mjj_" + to_string(fixed_cut_z) + "_ptmu_" + to_string(fixed_cut_w) + ".pdf");
  c_2D->Close();

// ---------------- Plotting 3D ----------------

  vector <TCanvas*> canvas_mu;
  for (Int_t m =0; m < bins[3]; m++)
    {
      canvas_mu.push_back(nullptr);
    }

  TString rootname = output+"/rate_3D_ptjet1_ptjet2_jetmass_ptmu.root";
  TFile *f=TFile::Open(rootname, "recreate");
  f->cd();

  for (Int_t m = 0; m < bins[3]; m++)
    {
      int num = xmin[3] + m;
      TString canvname = Form("pure_rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d", num);
      TString pngname = Form(output+"/pure_rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d.png", num);
      TString pdfname = Form(output+"/pure_rate_3D_ptjet1_ptjet2_jetmass_ptmu_%d.pdf", num);
      canvas_mu.at(m) = new TCanvas(canvname,canvname,900.,650.);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->SetTitle("");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->Draw("BOX2 Z");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitle("p_{T}^{Jet1}");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitleOffset(1.8);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitle("p_{T}^{Jet2}");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitleOffset(1.8);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitle("m_{jj}");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitleOffset(1.4);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMinimum(0);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMaximum(50);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->Write();
      gPad->Update();
      TPaletteAxis* pal = (TPaletteAxis*) ptjet1_ptjet2_jetmass_ptmu.at(m)->GetListOfFunctions()->FindObject("palette");
      pal->SetX1NDC(0.926); // It doesn't work (Why??)
      // ptjet1_ptjet2_jetmass_ptmu.at(m)->SaveAs(rootname);

      TLatex Tex4;
      Tex4.SetTextSize(0.035);
      Tex4.DrawLatexNDC(0.3,0.96,Form("Zero Bias L1 acceptance p_{T}^{#mu} > %d GeV", num));
      Tex4.Draw("same");

      canvas_mu.at(m)->SaveAs(pngname);
      canvas_mu.at(m)->SaveAs(pdfname);
      canvas_mu.at(m)->Close();
    }

  f->Close();

  rates_4D->SaveAs(output+"/PureRates_4D.root");

  cout << "\nNumber of bins = " << rates_4D->GetNbins() << endl;
  cout << "Number of good combinations = " << combinations << endl;

  Int_t Fixed_cut[4] = {1,1,1,1};
  cout << "\nRate for fixed cut {1,1,1,1} = " << rates_4D->GetBinContent(Fixed_cut) << "\n" << endl;

  return;
}
