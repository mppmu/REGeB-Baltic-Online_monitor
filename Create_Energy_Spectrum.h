//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  7 15:03:42 2020 by ROOT version 6.22/00
// from TTree data/Pulse Shape Output Tree
// found on file: htest.root
//////////////////////////////////////////////////////////

#ifndef Create_Energy_Spectrum_h
#define Create_Energy_Spectrum_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TLegend.h>

#include <iostream>
#include <fstream>


// Headers needed by this particular selector
#include <vector>
#include <string>

class SegBEGEK2ESA;

using std::vector;
using std::string;

class Create_Energy_Spectrum : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> info_idx = {fReader, "info_idx"};
   TTreeReaderValue<Double_t> info_time = {fReader, "info_time"};
   TTreeReaderArray<int> raw_pp_ch = {fReader, "raw_pp_ch"};
   TTreeReaderArray<int> raw_pp_trig_max = {fReader, "raw_pp_trig_max"};
   TTreeReaderArray<int> raw_pp_mca = {fReader, "raw_pp_mca"};
   TTreeReaderArray<int> raw_wf_smpl_v = {fReader, "raw_wf_smpl_v"};
   TTreeReaderArray<int> raw_wf_smpl_n = {fReader, "raw_wf_smpl_n"};


   Create_Energy_Spectrum(TTree * /*tree*/ =0) { }
   virtual ~Create_Energy_Spectrum() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Create_Energy_Spectrum,0);

private:
    TH1D* he_regeb;
    TH1D* he_xtrab;
    TH1F* hf_regeb;
    TH1F* hf_xtrab;
    float rough_regeb_energy_calibration_constant;
    float rough_xtrab_energy_calibration_constant;
    TString output1pngfilename;
    TString output2pngfilename;
    TString myoption;
   SegBEGeK2ESA* esa_tools;


};

#endif

#ifdef Create_Energy_Spectrum_cxx
void Create_Energy_Spectrum::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Create_Energy_Spectrum::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Create_Energy_Spectrum_cxx
