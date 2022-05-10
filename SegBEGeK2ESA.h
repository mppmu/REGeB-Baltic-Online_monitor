#ifndef __SegBEGeK2ESA__H
#define __SegBEGeK2ESA__H

/*
---------------------------------------------------------------------
  Class    : SegBEGeK2ESA
  Function : energy spectrum analysis
             take input energy spectrum, fit subrange with various functions

  Author   : Xiang Liu
  Date     : 01-12-2016
  Update     based on PXESA
---------------------------------------------------------------------
*/
/* C++ includes */
#include  <math.h>
#include  <iostream>
#include  <string.h>

/* ROOT MINUIT */
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TRandom2.h>
#include <TMinuit.h>
#include <TError.h>

/* ROOT includes */
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <string>
using std::string;

using namespace std;

const int maxnbins=50000;
const int maxnpeaks=20;

class SegBEGeK2ESA{

public:

// ==== cons and des ==== //
    SegBEGeK2ESA();
    SegBEGeK2ESA(TH1F* hmca);
    ~SegBEGeK2ESA();

// === setter and getter === //
    void SetInitialEscale(float escal) { escale_initial_guess=escal; };
    void SetInitialIntercept(float intercep) { intercept_initial_guess=intercep; };

//---> single peak fit
    float GetFittedSum()   { return fitted_sum;  };
    float GetFittedBkg()   { return fitted_bkg;  };
    float GetFittedMean()  { return fitted_mean; };
    float GetFittedSigma() { return fitted_sigma;};
    float GetFittedSumError()   { return fitted_sumerror;  };
    float GetFittedBkgError()   { return fitted_bkgerror;  };
    float GetFittedMeanError()  { return fitted_meanerror; };
    float GetFittedSigmaError() { return fitted_sigmaerror;};
//---> two gaussian peak fit
    float GetSecondFittedSum()   { return fitted_sum_two;  };
    float GetSecondFittedBkg()   { return fitted_bkg_two;  };
    float GetSecondFittedMean()  { return fitted_mean_two; };
    float GetSecondFittedSigma() { return fitted_sigma_two;};
    float GetSecondFittedSumError()   { return fitted_sumerror_two;  };
    float GetSecondFittedBkgError()   { return fitted_bkgerror_two;  };
    float GetSecondFittedMeanError()  { return fitted_meanerror_two; };
    float GetSecondFittedSigmaError() { return fitted_sigmaerror_two;};

//---> multipeak fit
    void SetNPeaks(int n) { npeaks=n; };
    void SetSpecPhotonEnergy(int n, float e) { fitspec_photone[n]=e; };
    void SetSpecPeakEnergy(int n, float e) { fitspec_peake[n]=e; };
    float GetEscale()         { return fittedspec_slope; };
    float GetEscaleError()    { return fittedspec_slope_error; };
    float GetIntercept()      { return fittedspec_intercept; };
    float GetInterceptError() { return fittedspec_intercept_error; };
    int   GetSpecNumPhotons()            { return npeaks; };
    float GetSpecPhotonEnergy(int n)     { return fitspec_photone[n]; };
    float GetSpecFittedMean(int n)       { return fittedspec_mean[n]; };
    float GetSpecFittedMeanError(int n)  { return fittedspec_meanerror[n]; };
    float GetSpecFittedBackground(int n)       { return fittedspec_bkg[n]; };
    float GetSpecFittedBackgroundError(int n)  { return fittedspec_bkgerror[n]; };
    float GetSpecFittedSum(int n)        { return fittedspec_sum[n]; };
    float GetSpecFittedSumError(int n)   { return fittedspec_sumerror[n]; };
    float GetSpecFittedSigma(int n)      { return fittedspec_sigma[n]; };
    float GetSpecFittedSigmaError(int n) { return fittedspec_sigmaerror[n]; };

// === fitting === //
    //---> this resolution is from REGe, too optimistic for other detectors
    float calculate_energy_resolution(float energy);
    int   get_index_from_histogram(TH1F* ht, float hval);
    void  set_histogram_to_be_fitted(TH1F* hmca, float escal, float intercep);
    int   get_index_after_setting_histogram(float hval);
    void  fit_subrange_gausspluspol(float e_min, float e_max, bool with_constraint);
    void  fit_subrange_twogausspluspol(float e_min, float e_max, float meanone, float meantwo);
    void  quality_check_fit_ba133_mca_spectrum(TCanvas* cfit);
    void  calibrate_ba133_mca_spectrum(TCanvas* cfit);
    void  calibrate_ba133_passivation_mca_spectrum(TCanvas* cfit);
    void  calibrate_mca_spectrum(TCanvas* cfit);
    void  calibrate_mca_spectrum_without_offset(TCanvas* cfit);

protected:

private:

//---> input original spectrum, take it as mca, it maybe calibrated, does not matter
   Int_t   mca_nbins;
   Float_t mca_entry[maxnbins]; // entry of each bin from the energy spectrum
   Float_t mca_binsize;
   Float_t mca_mine;
   Float_t mca_maxe;  // so mca_binsize=(mca_mine-mca_maxe)/mca_nbins

//---> take the initial guess on the energy scale and get the energy-related variables (unit keV)
//     fit is applied always to these variables
   Int_t   he_nbins;
   Float_t escale_initial_guess;
   Float_t intercept_initial_guess; // notice: energy=intercept+escale*mca
   Float_t he_entry[maxnbins];
   Float_t he_binsize;
   Float_t he_mine;
   Float_t he_maxe;
   TString he_name;

//---> fit a subrange of the input energy spectrum with gaussian-plus-1orderpolynomial function
//     Gaussian + a+b*x
   Float_t fit_peake;     // input variables
   Float_t fit_windowsize;
   Float_t fit_mine;
   Float_t fit_maxe;    
   Int_t   fit_minibin;
   Int_t   fit_maxibin;
   Int_t   fit_nbins;
   Float_t fitted_sum;       // output variables
   Float_t fitted_bkg;       // plus-minus 1.5*sigma
   Float_t fitted_mean;
   Float_t fitted_sigma;
   Float_t fitted_sumerror;
   Float_t fitted_bkgerror;
   Float_t fitted_meanerror;
   Float_t fitted_sigmaerror;
   Float_t fitted_a;
   Float_t fitted_b;
   Float_t fitted_aerror;
   Float_t fitted_berror;
   Float_t fitted_sum_two;  // these are when fitting with a double Gaussian function
   Float_t fitted_bkg_two;
   Float_t fitted_mean_two;
   Float_t fitted_sigma_two;
   Float_t fitted_sumerror_two;  
   Float_t fitted_bkgerror_two;
   Float_t fitted_meanerror_two;
   Float_t fitted_sigmaerror_two;

//---> list of energy peaks and fitted results, used when fitt a whole spectrum
//      and variables for calibration
   Int_t  npeaks;                    // input
   Float_t fitspec_photone[maxnpeaks]; // real-photon-energy
   Float_t fitspec_peake[maxnpeaks];   // peak-center-energy = real-photon-energy/energy_scale
   Float_t fitspec_windowsize[maxnpeaks];
   Float_t fitspec_mine[maxnpeaks];
   Float_t fitspec_maxe[maxnpeaks];
   Float_t fittedspec_sum[maxnpeaks];
   Float_t fittedspec_bkg[maxnpeaks];
   Float_t fittedspec_mean[maxnpeaks];
   Float_t fittedspec_sigma[maxnpeaks];
   Float_t fittedspec_sumerror[maxnpeaks];
   Float_t fittedspec_bkgerror[maxnpeaks];
   Float_t fittedspec_meanerror[maxnpeaks];
   Float_t fittedspec_sigmaerror[maxnpeaks];
//---> after fitting the fitspec_photone vs. fittedspec_mean, get energy scale etc
   Float_t fittedspec_slope;
   Float_t fittedspec_intercept;
   Float_t fittedspec_slope_error;
   Float_t fittedspec_intercept_error;


};

#endif
