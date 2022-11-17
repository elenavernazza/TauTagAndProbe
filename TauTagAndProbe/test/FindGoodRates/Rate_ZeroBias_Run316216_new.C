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
#include <string.h>
#include <vector>
#include <TPaletteAxis.h>
#include <TLatex.h>

using namespace std;

void Rate()
{

  Int_t bins[4] = {20, 20, 40, 13};
  Double_t xmin[4] = {0., 0., 200., 3.};
  Double_t xmax[4] = {100., 100., 1000., 16.};

  Int_t Denominator = 0;
  UInt_t noObjects_events = 0.;
  UInt_t GoodEvents = 0.;
  vector <UInt_t> GoodEvents_Mu (bins[3], 0);
  UInt_t noJets_events = 0.;
  UInt_t noLumi_events = 0.;
  UInt_t diffMass_events = 0.;

  cout << "Begin loop" << endl;

  TString path = "/grid_mnt/data__data.polcms/cms/vernazza/CMSSW_10_3_1/src/TauTagAndProbe/TauTagAndProbe/test/Run3_MC_VBFHToTauTau_M125_RAW/";

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

      inTree->SetBranchAddress("EventNumber", &in_EventNumber);
      inTree->SetBranchAddress("RunNumber", &in_RunNumber);
      inTree->SetBranchAddress("lumi", &in_lumi);
      inTree->SetBranchAddress("l1tPtJet", &in_l1tPtJet);
      inTree->SetBranchAddress("l1tEtaJet", &in_l1tEtaJet);
      inTree->SetBranchAddress("l1tPhiJet", &in_l1tPhiJet);
      inTree->SetBranchAddress("l1tMuPt", &in_l1tMuPt);
      inTree->SetBranchAddress("l1tMuEta", &in_l1tMuEta);
      inTree->SetBranchAddress("l1tMuPhi", &in_l1tMuPhi);

      cout << "Number of events in the Ntuple " << nt << " = " << inTree->GetEntries() << endl;

      for(UInt_t i = 0 ; i < inTree->GetEntries() ; ++i)
        {
          inTree->GetEntry(i);
          // if(i%10000==0) cout << "Entry #" << i << endl; 
          if(in_lumi<63 || in_lumi>250)
            {
              ++ noLumi_events;
              continue;
            }

          Float_t weight = 1.;

          ++Denominator;

          // Do a loop on the invariant mass for jets>20 and take the combination giving the highest mjj. The most energetic is jet1, the least energetic is the jet2

          bool atLeast2Jets = false;
          float highestMass = -1.;
          float L1Jet1Pt = -1.;
          float L1Jet2Pt = -1.;
          float Mass = -1.;

          if (in_l1tPtJet->size() > 0 && in_l1tMuPt->size() > 0)
            {

              for (UInt_t iL1Jet = 0 ; iL1Jet < in_l1tPtJet->size() ; ++iL1Jet)
                {
                  TLorentzVector L1Jet1;
                  L1Jet1.SetPtEtaPhiM(in_l1tPtJet->at(iL1Jet),in_l1tEtaJet->at(iL1Jet),in_l1tPhiJet->at(iL1Jet),0.);

                  for (UInt_t jL1Jet = 0 ; jL1Jet < in_l1tPtJet->size() ; ++jL1Jet)
                    {
                      TLorentzVector L1Jet2;
                      L1Jet2.SetPtEtaPhiM(in_l1tPtJet->at(jL1Jet),in_l1tEtaJet->at(jL1Jet),in_l1tPhiJet->at(jL1Jet),0.);
                      TLorentzVector DiJet = L1Jet1 + L1Jet2;
                      Mass = DiJet.M();

                      bool check = in_l1tPtJet->at(iL1Jet) > 20. && in_l1tPtJet->at(jL1Jet) > 20. && jL1Jet != iL1Jet ;
                      if (check)
                        {
                          atLeast2Jets = true;
                          if (Mass > highestMass) 
                            {
                              highestMass = Mass;
                              L1Jet1Pt = max(in_l1tPtJet->at(iL1Jet), in_l1tPtJet->at(jL1Jet));
                              L1Jet2Pt = min(in_l1tPtJet->at(iL1Jet), in_l1tPtJet->at(jL1Jet));
                            }
                        }
                    }  
                }       

          //     float subleadingPt = -1.;
          //     float Mass_1 = -1.;
          //     float highestJetPt = -1.;
          //     UInt_t highestJetPtIndex = -1; 

          //     for (UInt_t iL1Jet = 0 ; iL1Jet < in_l1tPtJet->size() ; ++iL1Jet)
          //       {
          //         if (in_l1tPtJet->at(iL1Jet) > highestJetPt) 
          //           {
          //             highestJetPt = in_l1tPtJet->at(iL1Jet);
          //             highestJetPtIndex = iL1Jet;
          //           }
          //       }

          // //     if (highestJetPtIndex != 0)
          // //       {
          // //         cout << "Highest Pt for event " << i << " is jet number " << highestJetPtIndex << " : " << in_l1tPtJet->at(highestJetPtIndex) << "  " << in_l1tPtJet->at(0) << endl;
          // //       }

          //     for (UInt_t iL1Jet = 0 ; iL1Jet < in_l1tPtJet->size() ; ++iL1Jet)
          //       {
          //         TLorentzVector L1Jet1;
          //         L1Jet1.SetPtEtaPhiM(in_l1tPtJet->at(iL1Jet),in_l1tEtaJet->at(iL1Jet),in_l1tPhiJet->at(iL1Jet),0.);

          //         for (UInt_t jL1Jet = 0 ; jL1Jet < in_l1tPtJet->size() && jL1Jet != iL1Jet ; ++jL1Jet)
          //           {
          //             TLorentzVector L1Jet2;
          //             L1Jet2.SetPtEtaPhiM(in_l1tPtJet->at(jL1Jet),in_l1tEtaJet->at(jL1Jet),in_l1tPhiJet->at(jL1Jet),0.);

          //             if(jL1Jet != highestJetPtIndex && L1Jet2.Pt()>subleadingPt) 
          //               {
          //                 subleadingPt = L1Jet2.Pt();
          //                 TLorentzVector DiJet = L1Jet1 + L1Jet2;
          //                 Mass_1 = DiJet.M();
          //               }

          //           }

          //       }

              float highestMuonPt = -1.;
              bool atLeast1Muon = false;

              for (UInt_t iL1Muon = 0 ; iL1Muon < in_l1tMuPt->size() ; ++iL1Muon)
                {
                  if (in_l1tMuPt->at(iL1Muon) > highestMuonPt) highestMuonPt = in_l1tMuPt->at(iL1Muon);
                  atLeast1Muon = true;	
                }

              if (atLeast2Jets && atLeast1Muon)
                {

                  PtJet1_PtJet2->Fill(L1Jet1Pt,L1Jet2Pt);
                  ++ GoodEvents;

                  // cout << "Event " << i << "  " << in_l1tPtJet->at(0) << "  " << in_l1tPtJet->at(1) << "  " << in_l1tMuPt->at(0) << endl;
                  // cout << "Event " << i << "  " << L1Jet1Pt << "  " << L1Jet2Pt << "  " << highestMuonPt << " " << highestMass << endl;
                  
                  if (L1Jet1Pt != in_l1tPtJet->at(0) || L1Jet2Pt != in_l1tPtJet->at(1))
                    {
                      TLorentzVector Jet1;
                      Jet1.SetPtEtaPhiM(in_l1tPtJet->at(0),in_l1tEtaJet->at(0),in_l1tPhiJet->at(0),0.);
                      TLorentzVector Jet2;
                      Jet2.SetPtEtaPhiM(in_l1tPtJet->at(1),in_l1tEtaJet->at(1),in_l1tPhiJet->at(1),0.);
                      TLorentzVector DiJets = Jet1 + Jet2;
                      float Mass_before = DiJets.M();
                      // cout << "Different: Mass before = " << Mass_before << " Mass after = " << highestMass << endl;
                      if (Mass_before < highestMass)
                        {
                         ++ diffMass_events;
                        }
                      else
                        {
                         cout << "Problem: Mass definition is wrong: Mass before = " << Mass_before << " Mass after = " << highestMass << endl;
                        }
                    }

                  if (L1Jet1Pt < L1Jet2Pt)
                    {
                      cout << "Problem: Jet1pt is smaller than Jet2pt: " << L1Jet1Pt << " " << L1Jet2Pt << endl;
                    }

                  for (Int_t n_cut = 0; n_cut < bins[3]; n_cut++)
                    {
                      int muon_cut_min = n_cut + xmin[3];
                      int muon_cut_max = n_cut + 1 + xmin[3];
                      if (n_cut == 0)
                        {
                          if (highestMuonPt < muon_cut_max)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(L1Jet1Pt, L1Jet2Pt, highestMass);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                      else if (n_cut == bins[3]-1)
                        {
                          if (highestMuonPt >= muon_cut_min)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(L1Jet1Pt, L1Jet2Pt, highestMass);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                      else
                        {
                          if (highestMuonPt >= muon_cut_min && highestMuonPt < muon_cut_max)
                            {
                              ptjet1_ptjet2_jetmass_ptmu.at(n_cut)->Fill(L1Jet1Pt, L1Jet2Pt, highestMass);
                              ++ GoodEvents_Mu.at(n_cut);
                            }
                        }
                    }
                }

              else
                {
                  // cout << "Event " << i << "  " << in_l1tPtJet->at(0) << "  " << in_l1tMuPt->at(0) << endl;
                  // cout << "Event " << i << "  Not checked\n" << endl;
                  ++ noJets_events;
                }

            }

          else
            {
              ++ noObjects_events;
              // cout << "Event " << i << "  No objects\n" << endl;
            }

        }

    }

  // for (int y=1; y<=bins[0]; y++)
  //   {
  //     for (int t=1; t<=bins[1]; t++)
  //       {
  //         for (int r=1; r<bins[2]; r++)
  //           {
  //             // Define a way to get the bin content in the 3 bins to check if we have events whe ptj1<ptj2
  //             if (y<t && ptjet1_ptjet2_jetmass_ptmu_15->GetBinContent())
  //               {
  //                 cout << "Jet1pt is smaller than Jet2pt: " << i_bin_x << " " << i_bin_y << " " << integral_xyzw << endl;
  //               }
  //           }
  //       }
  //   }

  float nb = 2544.;
  float thisLumiRun = 1.604E34;
  float scaleToLumi = 2.00E34;
  float scale = 0.001*(nb*11245.6)*scaleToLumi/thisLumiRun;

  cout << "Denominator = " << Denominator << endl;
  cout << "Number of events without jets or muons = " << noObjects_events << endl;
  cout << "Number of events without the second jet = " << noJets_events << endl;
  cout << "Number of events out of lumi = " << noLumi_events << endl;
  cout << "Number of events with different mass = " << diffMass_events << endl;
  cout << "Number of good events = " << GoodEvents << endl;

  for (Int_t i_n = 0; i_n < bins[3]; ++i_n)
    {
      cout << "Number of events " << i_n + xmin[3] << " = " << GoodEvents_Mu[i_n] << endl;
    }

  // Compute 4D histogram, where at each combination of jetPt1, jetPt2, Mjj and muonPt corresponds a rate

  THnF* rate_4D = new THnF("rates_4D","rates_4D", 4, bins, xmin, xmax);
  float integral_xyzw = -1.;
  int combinations = 0;

  for (Int_t i_bin_x = 1 ; i_bin_x <= bins[0] ; ++i_bin_x)
    {
      for (Int_t i_bin_y = 1 ; i_bin_y <= bins[1] ; ++i_bin_y)
        {
          for (Int_t i_bin_z = 1 ; i_bin_z <= bins[2] ; ++i_bin_z)
            {
              vector <float> int_xyz (bins[3], 0);
              // float int_xyz[bins[3]]
              for (Int_t w = 0; w < bins[3] ; w++)
                {
                  int_xyz.at(w) = ptjet1_ptjet2_jetmass_ptmu.at(w)->Integral(i_bin_x,bins[0]+1,i_bin_y,bins[1]+1,i_bin_z,bins[2]+1)/Denominator*scale;
                  // int_xyz[w] = int_ev;
                  // cout << "Integral scaled = " << int_ev << endl;
                }

              for (Int_t i_bin_w = 1 ; i_bin_w <= bins[3] ; ++i_bin_w)
                {
                  integral_xyzw = 0;
                  for (Int_t i_vec_w = 0; i_vec_w < bins[3]; ++i_vec_w)
                    {
                      if (i_vec_w >= i_bin_w - 1)
                        {
                          // integral_xyzw += int_xyz[i_vec_w];
                          integral_xyzw += int_xyz.at(i_vec_w);
                        }
                    }

                  if (integral_xyzw < 1.05 && integral_xyzw > 0.95 && i_bin_x >= i_bin_y)
                    {
                      // cout << "Bin " << i_bin_x << ", " << i_bin_y << ", " << i_bin_z << ", " << i_bin_w << " : " << integral_xyzw << endl;
                      ++ combinations;
                    }

                  if (i_bin_x == 1 && i_bin_y == 1 && i_bin_z == 1 && i_bin_w == 1)
                    {
                      cout << "Integral (1,1,1,1) = " << integral_xyzw << endl;
                    }
                  Int_t v_bin[4] = {i_bin_x, i_bin_y, i_bin_z, i_bin_w};
                  rates_4D->SetBinContent(v_bin, integral_xyzw);
                }
            }
        }
    }

  // --------------- Start ---------------
  // Compute 2D_rate with a fixed muonPt value and jetPt2 value 
  TH2F* rate_2D = new TH2F("rate_2D","rate_2D", bins[0], xmin[0], xmax[0], bins[1], xmin[1], xmax[1]);

  for (Int_t I_bin_x = 1 ; I_bin_x <= bins[0] ; ++I_bin_x)
    {
      for (Int_t I_bin_y = 1 ; I_bin_y <= bins[1] ; ++I_bin_y)
        {
          Float_t rate_value = 0;
          for (UInt_t i_mu = 6; i_mu < ptjet1_ptjet2_jetmass_ptmu.size(); ++i_mu)
            {
              rate_value += ptjet1_ptjet2_jetmass_ptmu.at(i_mu)->Integral(I_bin_x,bins[0]+1,I_bin_y,bins[1]+1,10,bins[2]+1)/Denominator*scale;
            }
          cout << rate_value << endl;
          rate_2D->SetBinContent(I_bin_x, I_bin_y, rate_value);
        }
    }

  TCanvas* c_2D = new TCanvas("c_2D","c_2D",700.,550.);
  c_2D->cd();
  gPad->SetPad(0.0,0.0,1.0,1.0);
  gPad->SetTopMargin(2);
  gPad->SetBottomMargin(2);
  gPad->SetRightMargin(2);
  rate_2D->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
  rate_2D->GetXaxis()->SetTitleOffset(1.3);
  // rate_2D->GetXaxis()->SetTitleSize(0.05);
  rate_2D->GetYaxis()->SetTitle("p_{T}^{j2} [GeV]");
  rate_2D->GetYaxis()->SetTitleOffset(1.3);
  rate_2D->GetZaxis()->SetTitle("rate [kHz]");
  rate_2D->GetZaxis()->SetTitleOffset(0.6);
  rate_2D->SetTitle("");
  // rate_2D->SetTitle("Zero Bias L1 rate M_{jj} > 420 GeV && p_{T}^{#mu} > 8 GeV");
  // rate_2D->SetTitleSize(0.3);
  // rate_2D->SetTitleOffset(2.8);
  rate_2D->SetStats(0);
  rate_2D->Draw("COLZ");

  TLatex Tex1;
  Tex1.SetTextSize(0.03);
  Tex1.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Simulation");
  Tex1.Draw("same");

  TLatex Tex2;
  Tex2.SetTextSize(0.035);
  Tex2.SetTextAlign(31);
  Tex2.DrawLatexNDC(0.90,0.91,"(14 TeV)");
  Tex2.Draw("same");

  TLatex Tex3;
  Tex3.SetTextSize(0.035);
  Tex3.DrawLatexNDC(0.25,0.96,"Zero Bias L1 rate M_{jj} > 420 GeV && p_{T}^{#mu} > 8 GeV");
  Tex3.Draw("same");

  c_2D->SaveAs("rate_ptj1_X_ptj2_Y_mjj_420_ptmu_8.png");
  c_2D->SaveAs("rate_ptj1_X_ptj2_Y_mjj_420_ptmu_8.pdf");
  c_2D->Close();

  // --------------- End ---------------

  for (int num = 0; num < bins[3]; num++)
    {
      Int_t Total_bin[4] = {1,1,1,num};
      cout << "Integral (1,1,1," << num << ") = " << rates_4D->GetBinContent(Total_bin) << endl;
    }

  cout << "Number of bins = " << rates_4D->GetNbins() << endl;
  cout << "Number of good combinations = " << combinations << endl;

  TCanvas* c = new TCanvas("c","c",800.,750.);

  PtJet1_PtJet2->Draw("COLZ");
  PtJet1_PtJet2->GetXaxis()->SetTitle("PtJet1");
  PtJet1_PtJet2->GetYaxis()->SetTitle("PtJet2");
  // for (int xx = 1; xx < bins[0]; ++xx)
  //   {
  //     for (int yy = 1; yy < bins[1]; ++yy)
  //       {
  //         int integrals = PtJet1_PtJet2->Integral(xx, bins[0], yy, bins[1]);
  //         if (xx < yy && integrals != 0)
  //           {
  //             cout << "Rate different from 0, binx " << xx << " biny " << yy << endl;
  //           }
  //       }
  //   }

  c->SaveAs("PtJet1_PtJet2.png");
  c->SaveAs("PtJet1_PtJet2.pdf");
  c->Close();

  vector <TCanvas*> canvas_mu;

  for (Int_t m =0; m < bins[3]; m++)
    {
      canvas_mu.push_back(nullptr);
    }

  for (Int_t m = 0; m < bins[3]; m++)
    {
      int num = xmin[3] + m;
      TString canvname = Form("ptjet1_ptjet2_jetmass_ptmu_%d", num);
      TString pngname = Form("ptjet1_ptjet2_jetmass_ptmu_%d.png", num);
      TString pdfname = Form("ptjet1_ptjet2_jetmass_ptmu_%d.pdf", num);
      TString rootname = Form("ptjet1_ptjet2_jetmass_ptmu_%d.root", num);
      canvas_mu.at(m) = new TCanvas(canvname,canvname,900.,650.);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->Draw("BOX2 Z");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitle("PtJet1");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetXaxis()->SetTitleOffset(1.8);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitle("PtJet2");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetYaxis()->SetTitleOffset(1.8);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitle("Mjj");
      ptjet1_ptjet2_jetmass_ptmu.at(m)->GetZaxis()->SetTitleOffset(1.4);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMinimum(0);
      ptjet1_ptjet2_jetmass_ptmu.at(m)->SetMaximum(200);
      gPad->Update();
      TPaletteAxis* pal = (TPaletteAxis*) ptjet1_ptjet2_jetmass_ptmu.at(m)->GetListOfFunctions()->FindObject("palette");
      pal->SetX1NDC(0.926); // It doesn't work (Why??)
      // ptjet1_ptjet2_jetmass_ptmu.at(m)->SaveAs(rootname);
      canvas_mu.at(m)->SaveAs(pngname);
      canvas_mu.at(m)->SaveAs(pdfname);
      canvas_mu.at(m)->Close();
    }

  rates_4D->SaveAs("Rates_4D.root");

  return;
}
