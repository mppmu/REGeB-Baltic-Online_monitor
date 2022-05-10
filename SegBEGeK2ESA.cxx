#include "SegBEGeK2ESA.h"

/* 
--------------------------------------------------
  following is function for fitting
  gaussian + polinormial 1
  ax+b+c*exp((x-d)^2/2e^2);
-------------------------------------------------- 
*/

Double_t mygaus(Double_t *x, Double_t *par){
  Double_t funcval=
         (par[0]/(2.0*3.14159*par[2]))
         *exp(0.0-(x[0]-par[1])*(x[0]-par[1])/(2.0*par[2]*par[2]));
  return funcval;
}


Double_t gauspluspol(Double_t *x, Double_t *par){
  Double_t funcval=par[0]+par[1]*x[0]
        +(par[2]/(2.0*3.14159*par[4]))
         *exp(0.0-(x[0]-par[3])*(x[0]-par[3])/(2.0*par[4]*par[4]));
  return funcval;
}

//---> y = b*x
Double_t mymca_vs_e(Double_t *x, Double_t *par)
{
   return par[0]*x[0];
}

/* 
--------------------------------------------------
  following is function for fitting
  gaussian1 + gaussian2 + polinormial 1
-------------------------------------------------- 
*/
Double_t twogauspluspol(Double_t *x, Double_t *par){
  Double_t funcval=par[0]+par[1]*x[0]
        +(par[2]/(sqrt(2.0*3.14159)*par[4]))
         *exp(0.0-(x[0]-par[3])*(x[0]-par[3])/(2.0*par[4]*par[4]))
        +(par[5]/(sqrt(2.0*3.14159)*par[7]))
         *exp(0.0-(x[0]-par[6])*(x[0]-par[6])/(2.0*par[7]*par[7]));
  return funcval;
}

// ---------------------------------------------------------------------------
SegBEGeK2ESA::SegBEGeK2ESA(){
   mca_nbins=0;
   he_nbins=0;
   escale_initial_guess=1.0;
   intercept_initial_guess=0.0;
   npeaks=0;
}

// ---------------------------------------------------------------------------
SegBEGeK2ESA::SegBEGeK2ESA(TH1F* hmca){
   mca_nbins=hmca->GetNbinsX();
   for (int i=0; i<mca_nbins; i++) mca_entry[i]=hmca->GetBinContent(i+1);
   mca_binsize=hmca->GetBinWidth(1);
   mca_mine   =hmca->GetBinCenter(1)        -mca_binsize*0.5;
   mca_maxe   =hmca->GetBinCenter(mca_nbins)+mca_binsize*0.5;
   escale_initial_guess=1.0;
   intercept_initial_guess=0.0;
   he_nbins=mca_nbins;
   for (int i=0; i<mca_nbins; i++) he_entry[i]=mca_entry[i];
   he_binsize = mca_binsize*escale_initial_guess;
   he_mine = mca_mine*escale_initial_guess+intercept_initial_guess;
   he_maxe = mca_maxe*escale_initial_guess+intercept_initial_guess;
   npeaks=0;
}

// ---------------------------------------------------------------------------
SegBEGeK2ESA::~SegBEGeK2ESA(){
;
}

// ---------------------------------------------------------------------------
float SegBEGeK2ESA::calculate_energy_resolution(float energy){
    Float_t eresolution =
    sqrt( 0.0405*0.0405*energy+
         (0.07726*energy/1000.0+1.084)*(0.07726*energy/1000.0+1.084) );
    return eresolution;
}

// ---------------------------------------------------------------------------
int SegBEGeK2ESA::get_index_from_histogram(TH1F* ht, Float_t hval)
{
    Int_t    nbins = ht->GetNbinsX();
    Float_t  binsize = ht->GetBinWidth(1);
    Float_t  fmin = ht->GetBinCenter(1)-0.5*binsize;
    //Float_t  fmax = ht->GetBinCenter(nbins)+0.5*binsize;
    Int_t    ih = int((hval-0.0001-fmin)/binsize)+1;
    //cout<<" energy "<<hval<<" ibin "<<ih<<endl;
    //cout<<fmin<<" "<<fmax<<" "<<binsize<<" "<<ih<<endl;
    return ih;
}

// ---------------------------------------------------------------------------
void SegBEGeK2ESA::set_histogram_to_be_fitted(TH1F* hmca, float escal, float intercep)
{
   mca_nbins=hmca->GetNbinsX();
   for (int i=0; i<mca_nbins; i++) mca_entry[i]=hmca->GetBinContent(i+1);
   mca_binsize=hmca->GetBinWidth(1);
   mca_mine   =hmca->GetBinCenter(1)        -mca_binsize*0.5;
   mca_maxe   =hmca->GetBinCenter(mca_nbins)+mca_binsize*0.5;
   escale_initial_guess=escal;
   intercept_initial_guess=intercep;
   he_nbins=mca_nbins;
   for (int i=0; i<mca_nbins; i++) he_entry[i]=mca_entry[i];
   he_binsize = mca_binsize*escale_initial_guess;
   he_mine = mca_mine*escale_initial_guess+intercept_initial_guess;
   he_maxe = mca_maxe*escale_initial_guess+intercept_initial_guess;
   he_name = TString(Form("%s",hmca->GetTitle()));
}

// ---------------------------------------------------------------------------
// this is the histogram bin index, so it starts from 1
// ---------------------------------------------------------------------------
int SegBEGeK2ESA::get_index_after_setting_histogram(Float_t hval)
{
     int ih=((hval-he_mine)/he_binsize)+1;
     if (ih<1) ih=1;
     if (ih>he_nbins) ih=he_nbins;
     return ih;
}

// ---------------------------------------------------------------------------
// here hmca is NOT necessarily calibrated, but e_min and e_max are real energy after calibration
// this is the reason why we first scale mca_?? to he_?? and then do the fit
// fit is always applied to he_?? variables
// this function must be called after set_histogram_to_be_fitted
// ---------------------------------------------------------------------------

void SegBEGeK2ESA::fit_subrange_gausspluspol(float e_min, float e_max, bool with_constraint)
{
   fit_peake = e_min+(e_max-e_min)*0.5;
   fitted_mean = fit_peake; // initial value
   fitted_sigma = calculate_energy_resolution(fit_peake); // initial value
   fit_windowsize = 3.0*fitted_sigma;
   fit_minibin = get_index_after_setting_histogram(e_min);
   fit_maxibin = get_index_after_setting_histogram(e_max);
   fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
   //fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize; // this is a bug
   fit_maxe = he_mine+float(fit_maxibin)*he_binsize;
   fit_nbins = fit_maxibin-fit_minibin+1;
   TH1F* histo_subrange = new TH1F("histo_subrange","subrange histogram",
                              fit_nbins, fit_mine, fit_maxe);
   fitted_sum=0.0;  // inivial value
   for (Int_t i=0; i<fit_nbins; i++) {
      histo_subrange->SetBinContent(i+1,he_entry[i+fit_minibin-1]);
      if (he_entry[i+fit_minibin-1]>fitted_sum) fitted_sum=he_entry[i+fit_minibin-1];
   }
   TF1* funcfit_singlepeak;
   funcfit_singlepeak = new TF1("funcfit_singlepeak",gauspluspol,fit_mine,fit_maxe,5);
   funcfit_singlepeak->SetParameter(0,histo_subrange->GetBinContent(1));
   funcfit_singlepeak->SetParameter(1,0.0);
   funcfit_singlepeak->SetParameter(2,histo_subrange->GetSum());
   funcfit_singlepeak->SetParameter(3,fit_peake);
   funcfit_singlepeak->SetParameter(4,fitted_sigma);
  if (with_constraint) {
   funcfit_singlepeak->SetParLimits(3,fit_peake-5.0,fit_peake+5.0);
   funcfit_singlepeak->SetParLimits(4,fitted_sigma*0.5,fitted_sigma*10.0);
  }

    histo_subrange->Fit(funcfit_singlepeak,"Q0NF");
    fitted_mean       = funcfit_singlepeak->GetParameter(3);
    fitted_meanerror  = funcfit_singlepeak->GetParError(3);
    fitted_sigma      = funcfit_singlepeak->GetParameter(4);
    fitted_sigmaerror = funcfit_singlepeak->GetParError(4);
   // fitted_sum        = funcfit_singlepeak->GetParameter(2);
   // fitted_sumerror   = funcfit_singlepeak->GetParError(2);
    fitted_a     = funcfit_singlepeak->GetParameter(0);
    fitted_aerror= funcfit_singlepeak->GetParError(0);
    fitted_b     = funcfit_singlepeak->GetParameter(1);
    fitted_berror= funcfit_singlepeak->GetParError(1);
    fitted_sum = 0.0;
    for (int i=1; i<=fit_nbins; i++)
      fitted_sum+=histo_subrange->GetBinContent(i)-
                  (fitted_a+fitted_b*histo_subrange->GetBinCenter(i));
    if (fitted_sum<0.0) fitted_sumerror=0.0;
    else                fitted_sumerror=sqrt(fitted_sum);
    int istart_bkg=int((fitted_mean-fitted_sigma*1.5-fit_mine)/he_binsize);
    int iend_bkg  =int((fitted_mean+fitted_sigma*1.5-fit_mine)/he_binsize);
    fitted_bkg = 0.0;
    for (int i=istart_bkg; i<=iend_bkg; i++) 
      fitted_bkg += histo_subrange->GetBinCenter(i)*fitted_b + fitted_a;
    if (fitted_bkg<0.0) fitted_bkgerror=0.0;
    else                fitted_bkgerror=sqrt(fitted_bkg);

    histo_subrange->SetLineColor(kBlack);
    histo_subrange->SetFillColor(kYellow);
    Float_t hismin=histo_subrange->GetMinimum();
    Float_t hismax=histo_subrange->GetMaximum();
    histo_subrange->SetMinimum(hismin*0.9);
    histo_subrange->SetMaximum(hismax*1.1);
    histo_subrange->SetTitle("");
    histo_subrange->SetStats(kFALSE);
    histo_subrange->SetXTitle("E[keV]");
    histo_subrange->SetYTitle("Entries");
    histo_subrange->SetTitleSize(0.07,"X");
    histo_subrange->SetTitleSize(0.07,"Y");
    histo_subrange->SetTitleOffset(0.7,"X");
    histo_subrange->SetTitleOffset(0.7,"Y");
    histo_subrange->SetLabelSize(0.05,"X");
    histo_subrange->SetLabelSize(0.04,"Y");
    histo_subrange->SetLabelOffset(0.001,"X");
    histo_subrange->SetLabelOffset(0.001,"Y");
    histo_subrange->DrawCopy();
    funcfit_singlepeak->SetLineColor(kRed);
    funcfit_singlepeak->SetLineWidth(1.0);
    funcfit_singlepeak->DrawCopy("same");

    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.08);
    l.DrawLatex(0.2,0.8,Form("mean %8.2f #pm %4.2f keV", fitted_mean, fitted_meanerror));
    l.DrawLatex(0.2,0.7,Form("FWHM %8.2f #pm %4.2f keV", fitted_sigma*2.35, fitted_sigmaerror*2.35));
    l.DrawLatex(0.2,0.6,Form("sum %8.2f", fitted_sum));

    histo_subrange->Delete();
    funcfit_singlepeak->Delete();
    return;
}

// ---------------------------------------------------------------------------
// this function must be called after set_histogram_to_be_fitted
// ---------------------------------------------------------------------------
void SegBEGeK2ESA::fit_subrange_twogausspluspol(float e_min, float e_max, float meanone, float meantwo)
{
   fit_peake = e_min+(e_max-e_min)*0.5;
   fitted_mean = meanone; // initial value
   fitted_mean_two = meantwo;
   fitted_sigma = calculate_energy_resolution(fit_peake); // initial value
   fitted_sigma_two = calculate_energy_resolution(fitted_mean_two);
   fit_windowsize = 3.0*fitted_sigma;
   fit_minibin = get_index_after_setting_histogram(e_min);
   fit_maxibin = get_index_after_setting_histogram(e_max);
   fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
   fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
   fit_nbins = fit_maxibin-fit_minibin+1;
   TH1F* histo_subrange = new TH1F("histo_subrange","subrange histogram",
                              fit_nbins, fit_mine, fit_maxe);
   fitted_sum=0.0;  // inivial value
   for (Int_t i=0; i<fit_nbins; i++) {
      histo_subrange->SetBinContent(i+1,he_entry[i+fit_minibin-1]);
      if (he_entry[i+fit_minibin-1]>fitted_sum) fitted_sum=he_entry[i+fit_minibin-1];
   }
   TF1* funcfit_twopeaks;
   funcfit_twopeaks = new TF1("funcfit_twopeaks",twogauspluspol,
                              fit_mine, fit_maxe, 8);

   funcfit_twopeaks->SetParameter(0,histo_subrange->GetBinContent(1));
   funcfit_twopeaks->SetParameter(1,0.0);
   funcfit_twopeaks->SetParameter(2,fitted_sum);
   funcfit_twopeaks->SetParLimits(2,0.0,1000000.0);
   funcfit_twopeaks->SetParameter(3,fitted_mean);
   funcfit_twopeaks->SetParLimits(3,fitted_mean-fitted_sigma,
                                      fitted_mean+fitted_sigma);
   funcfit_twopeaks->SetParameter(4,fitted_sigma);
   funcfit_twopeaks->SetParLimits(4,0.1*fitted_sigma,10.0*fitted_sigma);
   funcfit_twopeaks->SetParameter(5,fitted_sum);
   funcfit_twopeaks->SetParLimits(5,0.0,1000000.0);
   funcfit_twopeaks->SetParameter(6,fitted_mean_two);
   funcfit_twopeaks->SetParLimits(6,fitted_mean_two-fitted_sigma_two,
                                      fitted_mean_two+fitted_sigma_two);
   funcfit_twopeaks->SetParameter(7,fitted_sigma_two);
   funcfit_twopeaks->SetParLimits(7,0.1*fitted_sigma_two,10.0*fitted_sigma_two);

   histo_subrange->Fit(funcfit_twopeaks,"Q0N");
   fitted_mean      = funcfit_twopeaks->GetParameter(3);
   fitted_sigma     = funcfit_twopeaks->GetParameter(4);
   fitted_sum       = funcfit_twopeaks->GetParameter(2);
   fitted_meanerror = funcfit_twopeaks->GetParError(3);
   fitted_sigmaerror= funcfit_twopeaks->GetParError(4);
   fitted_a     = funcfit_twopeaks->GetParameter(0);
   fitted_b     = funcfit_twopeaks->GetParameter(1);
   fitted_aerror     = funcfit_twopeaks->GetParError(0);
   fitted_berror     = funcfit_twopeaks->GetParError(1);
   fitted_mean_two      = funcfit_twopeaks->GetParameter(6);
   fitted_sigma_two     = funcfit_twopeaks->GetParameter(7);
   fitted_sum_two       = funcfit_twopeaks->GetParameter(5);
   fitted_meanerror_two = funcfit_twopeaks->GetParError(6);
   fitted_sigmaerror_two= funcfit_twopeaks->GetParError(7);
//---> sum must be recalcualted
    fitted_sum = 0.0;
    for (int i=1; i<=fit_nbins; i++)
      fitted_sum+=histo_subrange->GetBinContent(i)-
                  (fitted_a+fitted_b*histo_subrange->GetBinCenter(i));
    if (fitted_sum<0.0) fitted_sumerror=0.0;
    else                fitted_sumerror=sqrt(fitted_sum);
    fitted_sum_two = fitted_sum;
    fitted_sumerror_two = fitted_sumerror;
//---> background
    int istart_bkg=int((fitted_mean-fitted_sigma*1.5-fit_mine)/he_binsize);
    int iend_bkg  =int((fitted_mean_two+fitted_sigma_two*1.5-fit_mine)/he_binsize);
    fitted_bkg = 0.0;
    for (int i=istart_bkg; i<=iend_bkg; i++) 
      fitted_bkg += histo_subrange->GetBinCenter(i)*fitted_b + fitted_a;
    if (fitted_bkg<0.0) fitted_bkgerror=0.0;
    else                fitted_bkgerror=sqrt(fitted_bkg);
    fitted_bkg_two = fitted_bkg;
    fitted_bkgerror_two = fitted_bkgerror;

   
    histo_subrange->SetLineColor(kBlack);
    histo_subrange->SetFillColor(kYellow);
    Float_t hismin=histo_subrange->GetMinimum();
    Float_t hismax=histo_subrange->GetMaximum();
    histo_subrange->SetMinimum(hismin*0.9);
    histo_subrange->SetMaximum(hismax*1.1);
    histo_subrange->SetTitle("");
    histo_subrange->SetStats(kFALSE);
    histo_subrange->SetXTitle("E(MeV)");
    histo_subrange->SetYTitle("Entries");
    histo_subrange->SetTitleSize(0.07,"X");
    histo_subrange->SetTitleSize(0.07,"Y");
    histo_subrange->SetTitleOffset(0.7,"X");
    histo_subrange->SetTitleOffset(0.7,"Y");
    histo_subrange->SetLabelSize(0.05,"X");
    histo_subrange->SetLabelSize(0.04,"Y");
    histo_subrange->SetLabelOffset(0.001,"X");
    histo_subrange->SetLabelOffset(0.001,"Y");
    histo_subrange->DrawCopy();
    funcfit_twopeaks->SetLineColor(kRed);
    funcfit_twopeaks->SetLineWidth(1.0);
    funcfit_twopeaks->DrawCopy("same");

    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.08);
    l.DrawLatex(0.2,0.8,Form("mean %8.4f pm %8.4f", fitted_mean, fitted_meanerror));
    l.DrawLatex(0.2,0.7,Form("FWHM %8.5f pm %8.5f", fitted_sigma*2.35, fitted_sigmaerror*2.35));
    l.DrawLatex(0.2,0.6,Form("sum %8.2f", fitted_sum));
    l.DrawLatex(0.2,0.5,Form("meantwo %8.4f pm %8.4f", fitted_mean_two, fitted_meanerror_two));
    l.DrawLatex(0.2,0.4,Form("FWHMtwo %8.5f pm %8.5f", fitted_sigma_two*2.35, fitted_sigmaerror_two*2.35));
    //l.DrawLatex(0.2,0.3,Form("sumtwo %8.2f", fitted_sum_two));

    histo_subrange->Delete();
    funcfit_twopeaks->Delete();
    return;
}

// ---------------------------------------------------------------------------
/*
 this function must be called after set_histogram_to_be_fitted
 this is a special version of fit_mca_spectrum_with_pdf
 the first photon line at 30 is in fact several gaussians with following mean values:

30.625      34.9 9      Cs Ka2
30.973      64.5 17     Cs Ka1
34.920      5.99 16     Cs Kb3
34.987      11.6 3      Cs Kb1
35.252      0.123 5     Cs Kb5
35.818      3.58 9      Cs Kb2
35.907      0.74 3      Cs Kb4
so the first peak is at (30.625*34.9+30.973*64.5)/(34.9+64.5)=30.851
the second peak ist at (34.920*5.99+34.987*11.6+35.252*0.123+35.818*3.58+35.907*0.74)/(5.99+11.6+0.123+3.58+0.74)=35.136

The following photon peaks are used
photon_energy_fwhm_fwhmerror_sum
30.851  1.0     1.0     1.0
35.136  1.0     1.0     1.0
80.9971 1.0     1.0     1.0
356.017 1.0     1.0     1.0

*/
// ---------------------------------------------------------------------------
void SegBEGeK2ESA::quality_check_fit_ba133_mca_spectrum(TCanvas* cfit)
{
    npeaks=4;
    fitspec_photone[0]=30.851;
    fitspec_photone[1]=35.136;
    fitspec_photone[2]=80.997;
    fitspec_photone[3]=356.017;

    for (int i=0; i<npeaks; i++) {
      fitspec_peake[i] = (fitspec_photone[i]-intercept_initial_guess)/escale_initial_guess;
      fitspec_windowsize[i] = 15.0/escale_initial_guess;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }

    cfit->cd();cfit->Divide(2,2);
    cfit->cd(1);
    fit_subrange_twogausspluspol(fitspec_mine[0], fitspec_maxe[1],
                                 fitspec_photone[0],fitspec_photone[1]);
    fittedspec_sum[0]        = fitted_sum;
    fittedspec_sumerror[0]   = fitted_sumerror;
    fittedspec_mean[0]       = fitted_mean;
    fittedspec_meanerror[0]  = fitted_meanerror;
    fittedspec_sigma[0]      = fitted_sigma;
    fittedspec_sigmaerror[0] = fitted_sigmaerror;
    fittedspec_sum[1]        = fitted_sum_two;
    fittedspec_sumerror[1]   = fitted_sumerror_two;
    fittedspec_mean[1]       = fitted_mean_two;
    fittedspec_meanerror[1]  = fitted_meanerror_two;
    fittedspec_sigma[1]      = fitted_sigma_two;
    fittedspec_sigmaerror[1] = fitted_sigmaerror_two;

    for (int ipeak=2; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(4); 
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(400.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->DrawCopy();
    cfit->Update(); 
    hpartenergy->Delete();

}

void SegBEGeK2ESA::calibrate_ba133_mca_spectrum(TCanvas* cfit)
{
    npeaks=7;
   fitspec_photone[0]=30.851;
   fitspec_photone[1]=35.136;
   fitspec_photone[2]=80.997;
   fitspec_photone[3]=276.398;
   fitspec_photone[4]=302.853;
   fitspec_photone[5]=356.017;
   fitspec_photone[6]=383.851;

    for (int i=0; i<npeaks; i++) {
      //fitspec_peake[i] = (fitspec_photone[i]-intercept_initial_guess)/escale_initial_guess;
      //fitspec_windowsize[i] = 15.0/escale_initial_guess;
      fitspec_peake[i] = fitspec_photone[i];
      fitspec_windowsize[i] = 14.0;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }

    cfit->cd();cfit->Divide(3,3);
    cfit->cd(1);
    fit_subrange_twogausspluspol(fitspec_mine[0], fitspec_maxe[1],
                                 fitspec_photone[0],fitspec_photone[1]);
    fittedspec_sum[0]        = fitted_sum;
    fittedspec_sumerror[0]   = fitted_sumerror;
    fittedspec_bkg[0]        = fitted_bkg;
    fittedspec_bkgerror[0]   = fitted_bkgerror;
    fittedspec_mean[0]       = fitted_mean;
    fittedspec_meanerror[0]  = fitted_meanerror;
    fittedspec_sigma[0]      = fitted_sigma;
    fittedspec_sigmaerror[0] = fitted_sigmaerror;
    fittedspec_sum[1]        = fitted_sum_two;
    fittedspec_sumerror[1]   = fitted_sumerror_two;
    fittedspec_bkg[1]        = fitted_bkg_two;
    fittedspec_bkgerror[1]   = fitted_bkgerror_two;
    fittedspec_mean[1]       = fitted_mean_two;
    fittedspec_meanerror[1]  = fitted_meanerror_two;
    fittedspec_sigma[1]      = fitted_sigma_two;
    fittedspec_sigmaerror[1] = fitted_sigmaerror_two;

    for (int ipeak=2; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      //fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kFALSE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_bkg[ipeak]        = fitted_bkg;
      fittedspec_bkgerror[ipeak]   = fitted_bkgerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(7); // linear fit of MCA vs energy
    float x[7];
    float ex[7];
    float y[7];
    float ey[7];
    for (int i=0; i<npeaks; i++) {
      x[i] = fittedspec_mean[i];
     ex[i] = fittedspec_meanerror[i];
      y[i] = fitspec_photone[i];
     ey[i] = 0.1;
    }
    TGraphErrors* gr_energy_vs_mca = new TGraphErrors(npeaks,x,y,ex,ey);
   //gr_energy_vs_mca->SetTitle("true-energy vs. pre-calibrated energy");
   gr_energy_vs_mca->SetTitle(Form("%s",he_name.Data()));
   gr_energy_vs_mca->SetMarkerStyle(21);
   gr_energy_vs_mca->GetXaxis()->SetTitle("MCA [ADC counts]");
   gr_energy_vs_mca->GetYaxis()->SetTitle("Energy [keV]");
   gr_energy_vs_mca->Draw("ALP");
/*
   gr_energy_vs_mca->Fit("pol1","Q0F");
   TF1* f_fit = gr_energy_vs_mca->GetFunction("pol1");
   fittedspec_slope = f_fit->GetParameter(1);
   fittedspec_intercept = f_fit->GetParameter(0);
   fittedspec_slope_error = f_fit->GetParError(1);
   fittedspec_intercept_error = f_fit->GetParError(0);
*/
   TF1* f_fit = new TF1("f_fit",mymca_vs_e,
                        x[0]-10.0,x[npeaks-1]+10.0, 1);
   gr_energy_vs_mca->Fit(f_fit,"Q0N");
   fittedspec_slope = f_fit->GetParameter(0);
   fittedspec_slope_error = f_fit->GetParError(0);
   fittedspec_intercept = 0.0;
   fittedspec_intercept_error = 0.0;
   f_fit->Draw("same");
   TLatex ll;
   ll.SetTextSize(0.08);
   ll.SetNDC();
   ll.DrawLatex(0.35,0.3,Form("y = %6.4f(%6.4f)",
                 fittedspec_intercept,fittedspec_intercept_error));
   ll.DrawLatex(0.35,0.2,Form(" + %6.5f(%6.5f) x",
                 fittedspec_slope,fittedspec_slope_error));

//--> residual
    cfit->cd(8);
    for (int i=0; i<npeaks; i++) {
       x[i] = fitspec_photone[i];
      ex[i] = 0.1;
       y[i] = fittedspec_mean[i]-fitspec_photone[i]/fittedspec_slope;
      ey[i] = fittedspec_meanerror[i];
    }
    TGraphErrors* gr_residual_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_residual_vs_energy->SetTitle("residual vs. true energy");
   gr_residual_vs_energy->SetMarkerStyle(21);
   gr_residual_vs_energy->GetXaxis()->SetTitle("E [keV]");
   gr_residual_vs_energy->GetYaxis()->SetTitle("residual [keV]");
   gr_residual_vs_energy->Draw("ALP");

//--> FWHM
    cfit->cd(9);
    for (int i=0; i<npeaks; i++) {
      x[i] = fitspec_photone[i];
     ex[i] = 0.1;
      y[i] = fittedspec_sigma[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
     ey[i] = fittedspec_sigmaerror[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
    }
    TGraphErrors* gr_fwhm_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
    gr_fwhm_vs_energy->SetTitle("resolution");
    gr_fwhm_vs_energy->SetMarkerStyle(21);
    gr_fwhm_vs_energy->GetXaxis()->SetTitle("Energy [keV]");
    gr_fwhm_vs_energy->GetYaxis()->SetTitle("FWHM [keV]");
    gr_fwhm_vs_energy->Draw("ALP");
/*
    cfit->cd(9);
    gPad->SetLogy();
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(400.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->GetXaxis()->SetTitle("E [keV]");
    hpartenergy->GetYaxis()->SetTitle("Entries");
    hpartenergy->Draw();
*/
    cfit->Update();

}

void SegBEGeK2ESA::calibrate_ba133_passivation_mca_spectrum(TCanvas* cfit)
{
    npeaks=4;
   fitspec_photone[0]=276.398;
   fitspec_photone[1]=302.853;
   fitspec_photone[2]=356.017;
   fitspec_photone[3]=383.851;

    for (int i=0; i<npeaks; i++) {
      fitspec_peake[i] = fitspec_photone[i];
      fitspec_windowsize[i] = 14.0;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }

    cfit->cd();cfit->Divide(3,3);

    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak+1);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      //fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kFALSE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_bkg[ipeak]        = fitted_bkg;
      fittedspec_bkgerror[ipeak]   = fitted_bkgerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(5); // linear fit of MCA vs energy
    float x[4];
    float ex[4];
    float y[4];
    float ey[4];
    for (int i=0; i<npeaks; i++) {
      x[i] = fittedspec_mean[i];
     ex[i] = fittedspec_meanerror[i];
      y[i] = fitspec_photone[i];
     ey[i] = 0.1;
    }
    TGraphErrors* gr_energy_vs_mca = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_energy_vs_mca->SetTitle("true-energy vs. pre-calibrated energy");
   gr_energy_vs_mca->SetMarkerStyle(21);
   gr_energy_vs_mca->GetXaxis()->SetTitle("MCA [ADC counts]");
   gr_energy_vs_mca->GetYaxis()->SetTitle("Energy [keV]");
   gr_energy_vs_mca->Draw("ALP");

   TF1* f_fit = new TF1("f_fit",mymca_vs_e,
                        x[0]-10.0,x[npeaks-1]+10.0, 1);
   gr_energy_vs_mca->Fit(f_fit,"Q0N");
   fittedspec_slope = f_fit->GetParameter(0);
   fittedspec_slope_error = f_fit->GetParError(0);
   fittedspec_intercept = 0.0;
   fittedspec_intercept_error = 0.0;
   f_fit->Draw("same");
   TLatex ll;
   ll.SetTextSize(0.08);
   ll.SetNDC();
   ll.DrawLatex(0.35,0.3,Form("y = %6.4f(%6.4f)",
                 fittedspec_intercept,fittedspec_intercept_error));
   ll.DrawLatex(0.35,0.2,Form(" + %6.5f(%6.5f) x",
                 fittedspec_slope,fittedspec_slope_error));

//--> residual
    cfit->cd(6);
    for (int i=0; i<npeaks; i++) {
       x[i] = fitspec_photone[i];
      ex[i] = 0.1;
       y[i] = fittedspec_mean[i]-fitspec_photone[i]/fittedspec_slope;
      ey[i] = fittedspec_meanerror[i];
    }
    TGraphErrors* gr_residual_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_residual_vs_energy->SetTitle("residual vs. true energy");
   gr_residual_vs_energy->SetMarkerStyle(21);
   gr_residual_vs_energy->GetXaxis()->SetTitle("E [keV]");
   gr_residual_vs_energy->GetYaxis()->SetTitle("residual [keV]");
   gr_residual_vs_energy->Draw("ALP");

//--> FWHM
    cfit->cd(7);
    for (int i=0; i<npeaks; i++) {
      x[i] = fitspec_photone[i];
     ex[i] = 0.1;
      y[i] = fittedspec_sigma[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
     ey[i] = fittedspec_sigmaerror[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
    }
    TGraphErrors* gr_fwhm_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
    gr_fwhm_vs_energy->SetTitle("resolution");
    gr_fwhm_vs_energy->SetMarkerStyle(21);
    gr_fwhm_vs_energy->GetXaxis()->SetTitle("Energy [keV]");
    gr_fwhm_vs_energy->GetYaxis()->SetTitle("FWHM [keV]");
    gr_fwhm_vs_energy->Draw("ALP");

//--> energy spectrum from 0 to 400
    cfit->cd(8);
    gPad->SetLogy();
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(400.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->GetXaxis()->SetTitle("E [keV]");
    hpartenergy->GetYaxis()->SetTitle("Entries");
    hpartenergy->Draw();


    cfit->Update();

}


/*
  this function must be called after setting number of photon lines and photon energies 
*/
void SegBEGeK2ESA::calibrate_mca_spectrum(TCanvas* cfit)
{

    for (int i=0; i<npeaks; i++) {
      fitspec_peake[i] = fitspec_photone[i]; // already considered energy-scale
      fitspec_windowsize[i] = 15.0;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }

    cfit->cd();cfit->Divide(2,4);
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak+1);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_bkg[ipeak]        = fitted_bkg;
      fittedspec_bkgerror[ipeak]   = fitted_bkgerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(npeaks+1); // linear fit of MCA vs energy
    float x[5];
    float ex[5];
    float y[5];
    float ey[5];
    for (int i=0; i<npeaks; i++) {
      x[i] = fittedspec_mean[i];
     ex[i] = fittedspec_meanerror[i];
      y[i] = fitspec_photone[i];
     ey[i] = 0.1;
    }
    TGraphErrors* gr_energy_vs_mca = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_energy_vs_mca->SetTitle("true-energy vs. pre-calibrated energy");
   gr_energy_vs_mca->SetMarkerStyle(21);
   gr_energy_vs_mca->GetXaxis()->SetTitle("MCA [ADC counts]");
   gr_energy_vs_mca->GetYaxis()->SetTitle("Energy [keV]");
   gr_energy_vs_mca->Draw("ALP");
   gr_energy_vs_mca->Fit("pol1","Q0F");
   TF1* f_fit = gr_energy_vs_mca->GetFunction("pol1");
   fittedspec_slope = f_fit->GetParameter(1);
   fittedspec_intercept = f_fit->GetParameter(0);
   fittedspec_slope_error = f_fit->GetParError(1);
   fittedspec_intercept_error = f_fit->GetParError(0);
   f_fit->Draw("same");
   TLatex ll;
   ll.SetTextSize(0.08);
   ll.SetNDC();
   ll.DrawLatex(0.35,0.3,Form("y = %6.4f(%6.4f)",
                 fittedspec_intercept,fittedspec_intercept_error));
   ll.DrawLatex(0.35,0.2,Form(" + %6.5f(%6.5f) x",
                 fittedspec_slope,fittedspec_slope_error));

    cfit->cd(npeaks+2);
    for (int i=0; i<npeaks; i++) {
      x[i] = fitspec_photone[i];
     ex[i] = 0.1;
      y[i] = fittedspec_sigma[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
     ey[i] = fittedspec_sigmaerror[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
    }
    TGraphErrors* gr_fwhm_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
    gr_fwhm_vs_energy->SetTitle("resolution");
    gr_fwhm_vs_energy->SetMarkerStyle(21);
    gr_fwhm_vs_energy->GetXaxis()->SetTitle("Energy [keV]");
    gr_fwhm_vs_energy->GetYaxis()->SetTitle("FWHM [keV]");
    gr_fwhm_vs_energy->Draw("ALP");

    cfit->cd(npeaks+3);
    gPad->SetLogy();
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(3000.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->GetXaxis()->SetTitle("E [keV]");
    hpartenergy->GetYaxis()->SetTitle("Entries");
    hpartenergy->Draw();
    cfit->Update();

}
/*
void SegBEGeK2ESA::calibrate_mca_spectrum_without_offset(TCanvas* cfit)
{

    for (int i=0; i<npeaks; i++) {
      fitspec_peake[i] = fitspec_photone[i]; // already considered energy-scale
      fitspec_windowsize[i] = 15.0;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }

    cfit->cd();cfit->Divide(2,3);
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak+1);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_bkg[ipeak]        = fitted_bkg;
      fittedspec_bkgerror[ipeak]   = fitted_bkgerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(npeaks+1); // linear fit of MCA vs energy
    float x[10];
    float ex[10];
    float y[10];
    float ey[10];
    for (int i=0; i<npeaks; i++) {
      x[i] = fittedspec_mean[i];
     ex[i] = fittedspec_meanerror[i];
      y[i] = fitspec_photone[i];
     ey[i] = 0.1;
    }
    TGraphErrors* gr_energy_vs_mca = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_energy_vs_mca->SetTitle("true-energy vs. pre-calibrated energy");
   gr_energy_vs_mca->SetMarkerStyle(21);
   gr_energy_vs_mca->GetXaxis()->SetTitle("MCA [ADC counts]");
   gr_energy_vs_mca->GetYaxis()->SetTitle("Energy [keV]");
   gr_energy_vs_mca->Draw("ALP");

   TF1* funcfit_myfit;
   funcfit_myfit = new TF1("funcfit_myfit",mymca_vs_e,
                              x[0]-10.0,x[npeaks-1]+10.0, 1);
   funcfit_myfit->SetParameter(0,1.0);

   //gr_energy_vs_mca->Fit("pol1","Q0F");
   //TF1* f_fit = gr_energy_vs_mca->GetFunction("pol1");
   //fittedspec_slope = f_fit->GetParameter(1);
   //fittedspec_intercept = f_fit->GetParameter(0);
   //fittedspec_slope_error = f_fit->GetParError(1);
   //fittedspec_intercept_error = f_fit->GetParError(0);

   gr_energy_vs_mca->Fit(funcfit_myfit,"Q0N");
   fittedspec_slope = funcfit_myfit->GetParameter(0);
   fittedspec_slope_error = funcfit_myfit->GetParError(0);
   fittedspec_intercept = 0.0;
   fittedspec_intercept_error = 0.0;
   funcfit_myfit->Draw("same");
   TLatex ll;
   ll.SetTextSize(0.08);
   ll.SetNDC();
   ll.DrawLatex(0.35,0.3,Form("y = %6.4f(%6.4f)",
                 fittedspec_intercept,fittedspec_intercept_error));
   ll.DrawLatex(0.35,0.2,Form(" + %6.5f(%6.5f) x",
                 fittedspec_slope,fittedspec_slope_error));

//--> residual
    cfit->cd(npeaks+2);
    for (int i=0; i<npeaks; i++) {
       x[i] = fitspec_photone[i];
      ex[i] = 0.1;
       y[i] = fittedspec_mean[i]-fitspec_photone[i]/fittedspec_slope;
      ey[i] = fittedspec_meanerror[i];
    }
    TGraphErrors* gr_residual_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_residual_vs_energy->SetTitle("residual vs. true energy");
   gr_residual_vs_energy->SetMarkerStyle(21);
   gr_residual_vs_energy->GetXaxis()->SetTitle("E [keV]");
   gr_residual_vs_energy->GetYaxis()->SetTitle("residual [keV]");
   gr_residual_vs_energy->Draw("ALP");

//---> FWHM
    cfit->cd(npeaks+3);
    for (int i=0; i<npeaks; i++) {
      x[i] = fitspec_photone[i];
     ex[i] = 0.1;
      y[i] = fittedspec_sigma[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
     ey[i] = fittedspec_sigmaerror[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
    }
    TGraphErrors* gr_fwhm_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
    gr_fwhm_vs_energy->SetTitle("resolution");
    gr_fwhm_vs_energy->SetMarkerStyle(21);
    gr_fwhm_vs_energy->GetXaxis()->SetTitle("Energy [keV]");
    gr_fwhm_vs_energy->GetYaxis()->SetTitle("FWHM [keV]");
    gr_fwhm_vs_energy->Draw("ALP");

//---> calibrated energy spectrum
    cfit->cd(npeaks+4);
    gPad->SetLogy();
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(3000.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->GetXaxis()->SetTitle("E [keV]");
    hpartenergy->GetYaxis()->SetTitle("Entries");
    hpartenergy->Draw();
    cfit->Update();

}
*/
void SegBEGeK2ESA::calibrate_mca_spectrum_without_offset(TCanvas* cfit)
{

    for (int i=0; i<npeaks; i++) {
      fitspec_peake[i] = fitspec_photone[i]; // already considered energy-scale
      fitspec_windowsize[i] = 15.0;
      fitspec_mine[i]=fitspec_peake[i]-fitspec_windowsize[i];
      fitspec_maxe[i]=fitspec_peake[i]+fitspec_windowsize[i];
    }
        int h_columns = 2;//<- modification by martin, so there are always enough subcanvases for given number of photopeaks plus the 4 other plots.
	int h_rows= (npeaks+4)/h_columns;// -"-
	if ((npeaks+4)%h_columns >0){h_rows++;}// -"-
    cfit->cd();cfit->Divide(h_columns,h_rows);// -"-
    
    for (int ipeak=0; ipeak<npeaks; ipeak++) {
      cfit->cd(ipeak+1);
      fit_subrange_gausspluspol(fitspec_mine[ipeak],fitspec_maxe[ipeak],kTRUE);
      fittedspec_sum[ipeak]        = fitted_sum;
      fittedspec_sumerror[ipeak]   = fitted_sumerror;
      fittedspec_bkg[ipeak]        = fitted_bkg;
      fittedspec_bkgerror[ipeak]   = fitted_bkgerror;
      fittedspec_mean[ipeak]       = fitted_mean;
      fittedspec_meanerror[ipeak]  = fitted_meanerror;
      fittedspec_sigma[ipeak]      = fitted_sigma;
      fittedspec_sigmaerror[ipeak] = fitted_sigmaerror;
    }

    cfit->cd(npeaks+1); // linear fit of MCA vs energy
    float x[10];
    float ex[10];
    float y[10];
    float ey[10];
    for (int i=0; i<npeaks; i++) {
      x[i] = fittedspec_mean[i];
     ex[i] = fittedspec_meanerror[i];
      y[i] = fitspec_photone[i];
     ey[i] = 0.1;
    }
    TGraphErrors* gr_energy_vs_mca = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_energy_vs_mca->SetTitle("true-energy vs. pre-calibrated energy");
   gr_energy_vs_mca->SetMarkerStyle(21);
   gr_energy_vs_mca->GetXaxis()->SetTitle("MCA [ADC counts]");
   gr_energy_vs_mca->GetYaxis()->SetTitle("Energy [keV]");
   gr_energy_vs_mca->Draw("ALP");

   TF1* funcfit_myfit;
   funcfit_myfit = new TF1("funcfit_myfit",mymca_vs_e,
                              x[0]-10.0,x[npeaks-1]+10.0, 1);
   funcfit_myfit->SetParameter(0,1.0);
/*
   gr_energy_vs_mca->Fit("pol1","Q0F");
   TF1* f_fit = gr_energy_vs_mca->GetFunction("pol1");
   fittedspec_slope = f_fit->GetParameter(1);
   fittedspec_intercept = f_fit->GetParameter(0);
   fittedspec_slope_error = f_fit->GetParError(1);
   fittedspec_intercept_error = f_fit->GetParError(0);
*/
   gr_energy_vs_mca->Fit(funcfit_myfit,"Q0N");
   fittedspec_slope = funcfit_myfit->GetParameter(0);
   fittedspec_slope_error = funcfit_myfit->GetParError(0);
   fittedspec_intercept = 0.0;
   fittedspec_intercept_error = 0.0;
   funcfit_myfit->Draw("same");
   TLatex ll;
   ll.SetTextSize(0.08);
   ll.SetNDC();
   ll.DrawLatex(0.35,0.3,Form("y = %6.4f(%6.4f)",
                 fittedspec_intercept,fittedspec_intercept_error));
   ll.DrawLatex(0.35,0.2,Form(" + %6.5f(%6.5f) x",
                 fittedspec_slope,fittedspec_slope_error));

//--> residual
    cfit->cd(npeaks+2);
    for (int i=0; i<npeaks; i++) {
       x[i] = fitspec_photone[i];
      ex[i] = 0.1;
       y[i] = fittedspec_mean[i]-fitspec_photone[i]/fittedspec_slope;
      ey[i] = fittedspec_meanerror[i];
    }
    TGraphErrors* gr_residual_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
   gr_residual_vs_energy->SetTitle("residual vs. true energy");
   gr_residual_vs_energy->SetMarkerStyle(21);
   gr_residual_vs_energy->GetXaxis()->SetTitle("E [keV]");
   gr_residual_vs_energy->GetYaxis()->SetTitle("residual [keV]");
   gr_residual_vs_energy->Draw("ALP");

//---> FWHM
    cfit->cd(npeaks+3);
    for (int i=0; i<npeaks; i++) {
      x[i] = fitspec_photone[i];
     ex[i] = 0.1;
      y[i] = fittedspec_sigma[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
     ey[i] = fittedspec_sigmaerror[i]*2.35/fittedspec_mean[i]*fitspec_photone[i];
    }
    TGraphErrors* gr_fwhm_vs_energy = new TGraphErrors(npeaks,x,y,ex,ey);
    gr_fwhm_vs_energy->SetTitle("resolution");
    gr_fwhm_vs_energy->SetMarkerStyle(21);
    gr_fwhm_vs_energy->GetXaxis()->SetTitle("Energy [keV]");
    gr_fwhm_vs_energy->GetYaxis()->SetTitle("FWHM [keV]");
    gr_fwhm_vs_energy->Draw("ALP");

//---> calibrated energy spectrum
    cfit->cd(npeaks+4);
    gPad->SetLogy();
    fit_minibin = get_index_after_setting_histogram(0.0);
    fit_maxibin = get_index_after_setting_histogram(3000.0);
    fit_nbins   = fit_maxibin-fit_minibin+1;
    fit_mine = he_mine+float(fit_minibin-1)*he_binsize;
    fit_maxe = he_mine+float(fit_maxibin-1)*he_binsize;
    TH1F* hpartenergy = new TH1F("hpartenergy","hpartenergy",fit_nbins, fit_mine, fit_maxe);
    for (int i=0; i<fit_nbins; i++) hpartenergy->SetBinContent(i+1,he_entry[i+fit_minibin]);
    hpartenergy->SetTitle(he_name.Data());
    hpartenergy->GetXaxis()->SetTitle("E [keV]");
    hpartenergy->GetYaxis()->SetTitle("Entries");
    hpartenergy->Draw();
    cfit->Update();

}
