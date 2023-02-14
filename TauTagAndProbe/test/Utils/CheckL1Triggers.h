#ifndef CHECKL1RIGGERS_H
#define CHECKL1RIGGERS_H

bool CheckGoodMuon (bool in_MuIso, int in_MuID);

bool CheckGoodJet (int in_JetID, TString JetIDType);

void CheckVBF (TTree* inTree, 
                UInt_t  i_ev, 
                vector<array<Float_t, 4>> set_of_on_cuts, 
                vector<array<Float_t, 4>> set_of_off_cuts, 
                bool pass_MuTau, vector<UInt_t>* acceptance_VBF, 
                vector<UInt_t>* acceptance_MuTau_VBF, 
                TString JetIDType, 
                TString Method, 
                bool JetSel30);

bool CheckMuTau (TTree* inTree, 
                UInt_t i_ev, 
                TString JetIDType, 
                TString Method, 
                bool JetSel30);

void FindGoodRates(THnF* rate_4D, 
                vector<array<Float_t, 4>>* set_of_on_cuts, 
                vector<array<Float_t, 4>>* set_of_off_cuts, 
                vector<array<Int_t, 4>>* set_of_on_bins, 
                Int_t starting_x_bin);

void PrintMuons (TTree* inTree, UInt_t  i_ev);

void PlotGain_2D_ptj1_mjj (Float_t fixed_cut_y, 
                Float_t fixed_cut_w, 
                Int_t bin0, Int_t xmin0, Int_t xmax0, 
                Int_t bin2, Int_t xmin2, Int_t xmax2, 
                vector<array<Float_t, 4>> set_of_on_cuts, 
                vector<array<Int_t, 4>> set_of_on_bins, 
                vector<UInt_t> acceptance_VBF, 
                vector<UInt_t> acceptance_MuTau_VBF, 
                UInt_t acceptance_MuTau, 
                TString Output);

#endif