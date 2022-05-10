#define Create_Energy_Spectrum_cxx
// The class definition in Create_Energy_Spectrum.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Create_Energy_Spectrum.C")
// root> T->Process("Create_Energy_Spectrum.C","some options")
// root> T->Process("Create_Energy_Spectrum.C+")
//

#include "SegBEGeK2ESA.h"
#include "Create_Energy_Spectrum.h"
#include <TH2.h>
#include <TStyle.h>

void Create_Energy_Spectrum::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   //rough_xtrab_energy_calibration_constant=1460.5/1.98805e+06; // starting 28.10.2020, for GCDX, 1500-risetime 550-flattop
   //rough_xtrab_energy_calibration_constant=1460.5/2.64875e+06; // starting 30.03.2021, for GCDX, 2000-risetime 250-flattop
   //rough_xtrab_energy_calibration_constant=2614.5/2.37382e+06; // starting 12.05.2021, for GCDX, 1000-risetime 200-flattop
   //rough_xtrab_energy_calibration_constant=1460.5/5.23046e+05; // starting 18.01.2022, for GCDX, 1000-risetime 200-flattop
   rough_xtrab_energy_calibration_constant=1460.5/1.31972e+06; // starting 04.04.2022, for GCDX, 1000-risetime 200-flattop
   rough_regeb_energy_calibration_constant=1460.5/4.46767e+06; // starting 28.10.2020, for REGeB
   //rough_xtrab_energy_calibration_constant=1460.5/2.00707e+06; // starting 07.10.2020, for GCDX, 1500-risetime 550-flattop
   //rough_regeb_energy_calibration_constant=1460.5/4.49992e+06; // starting 07.10.2020, for REGeB
   //rough_xtrab_energy_calibration_constant=1460.5/2.67652e+06; // starting 07.10.2020, for GCDX
   //rough_regeb_energy_calibration_constant=1460.5/5.98494e+06; // starting 07.10.2020, for REGeB
   he_regeb = new TH1D("REGeB_Energy","REGeB_Energy",
                       3100,-100.0,3000.0);
   he_xtrab = new TH1D("XtraB_Energy","XtraB_Energy",
                       3100,-100.0,3000.0);
   hf_regeb = new TH1F("fREGeB_Energy","fREGeB_Energy",
                       3100,-100.0,3000.0);
   hf_xtrab = new TH1F("fXtraB_Energy","fXtraB_Energy",
                       3100,-100.0,3000.0);
   output1pngfilename=TString("plot_regeb_latest.png");
   output2pngfilename=TString("plot_regeb_")
                     +option
                     +TString(".png");
   myoption=option;
   esa_tools = new SegBEGeK2ESA();
}

void Create_Energy_Spectrum::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Create_Energy_Spectrum::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   //cout<<" now event "<<entry<<" "<<endl;

   fReader.SetLocalEntry(entry);

//since 2022.01.18, put Baltic to channel 1 on minidex-fadc
   if (raw_pp_ch[0]==1) {
      he_xtrab->Fill( float(raw_pp_mca[0]) * rough_xtrab_energy_calibration_constant );
      hf_xtrab->Fill( float(raw_pp_mca[0]) * rough_xtrab_energy_calibration_constant );
   }
   else if (raw_pp_ch[0]==2) {
      he_regeb->Fill( float(raw_pp_mca[0]) * rough_regeb_energy_calibration_constant );
      hf_regeb->Fill( float(raw_pp_mca[0]) * rough_regeb_energy_calibration_constant );
   }


   return kTRUE;
}

void Create_Energy_Spectrum::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Create_Energy_Spectrum::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
//
//---> open the 2017 nov background only
//
   TFile* flastnov = new TFile("sum_histograms/h.regeb.2017nov.root");
   TH1F* hnov = (TH1F*) gDirectory->Get("he");
   hnov->Rebin(10);
   hnov->Scale(1.0/404.0);
   TLatex lt;
   lt.SetNDC();
   lt.SetTextSize(0.10);
//
   TCanvas* ce = new TCanvas("ce","ce",10,10,900,1200);
   ce->cd();
   ce->Divide(1,3);
   ce->cd(1);
   gPad->SetLogy();
   float ymax_nov   =hnov    ->GetMaximum();
   float ymax_xtrab =he_xtrab->GetMaximum();
   float ymax_regeb =he_regeb->GetMaximum();
   if (ymax_xtrab<ymax_regeb) ymax_xtrab=ymax_regeb;
   if (ymax_xtrab<ymax_nov  ) ymax_xtrab=ymax_nov;
   he_xtrab->SetMaximum(2.0*ymax_xtrab);
   he_xtrab->SetMinimum(0.5);
   he_xtrab->SetTitle(Form("Energy Spectrum %s",myoption.Data()));
   he_xtrab->GetXaxis()->SetTitle("Energy [keV]");
   he_xtrab->GetYaxis()->SetTitle("Events per hour keV");
   he_xtrab->SetStats(kFALSE);
   he_xtrab->SetLineColor(kBlue);
   he_xtrab->Draw();
   he_regeb->SetLineColor(kBlack);
   he_regeb->Draw("same");
   hnov->SetLineColor(kRed);
   hnov->Draw("same");
   TLegend* ll = new TLegend(0.45,0.6,0.88,0.88);
   ll->SetFillStyle(0);
   ll->AddEntry(he_xtrab,"GCDX");
   ll->AddEntry(he_regeb,"REGeB");
   ll->AddEntry(hnov,"REGeB Nov. 2017");
   ll->Draw();
   ce->cd(2);
   esa_tools->set_histogram_to_be_fitted(hf_xtrab,1.0,0.0);
   esa_tools->fit_subrange_gausspluspol(1440.0,1480.0,kFALSE);
   lt.DrawLatex(0.7,0.7,"GCDX");
   ce->cd(3);
   esa_tools->set_histogram_to_be_fitted(hf_regeb,1.0,0.0);
   esa_tools->fit_subrange_gausspluspol(1440.0,1480.0,kFALSE);
   lt.DrawLatex(0.7,0.7,"REGeB");
   ce->Update();
   ce->Print(output1pngfilename.Data());
   ce->Print(output2pngfilename.Data());

}
