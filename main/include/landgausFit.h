#pragma once
#include "Rtypes.h"


// incomplete function types to prevent pollution of the Header file


#define ERROR_RETURN_VALUE -1
#define SUCCESS_RETURN_VALUE 1

#ifndef __CINT__
#define DLL_exp _declspec(dllexport) 

class TF1;
class TH1D;
class TH2D;
class TGraph;
class string;

#else 
// if one rebuilds the dict files one has to copy this includes to the dict.cxx files otherwise it will not work 

#include <fstream>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#define DLL_exp 
#endif // __CINT__

// TGraph DLL_exp Analyse(TH2D* effi, Int_t size);
// TH2D DLL_exp loopOverFolder( Int_t &size);
// TGraph DLL_exp getStrip(TH2D *effi,Int_t stripNumber);
// TH2D DLL_exp load_ascii3(const char* filename);
// void DLL_exp save2ascii(const char* name, const TH2D& inputHist );

class DLL_exp landgausFit
{
public:
	landgausFit(void);
	~landgausFit(void);

	//////////////////////////////////////////////////////////////////////////
	// Setters 
	void setStartValues(Double_t Amplitude=1, Double_t landau_mean=200,Double_t landau_sigma=70,Double_t gauss_sigma=10);
	void setStartAmplitude(Double_t Amplitude=1);
	void setStartLandauMean(Double_t landau_mean=200);
	void setStartLandauSigma(Double_t landau_sigma=70);
	void setStartGaussSigma(Double_t gaus_sigma=10);
	void setLowerLimits(Double_t Amplitude=0, Double_t landau_mean=0,Double_t landau_sigma=0,Double_t gauss_sigma=0);
	void setUpperLimits(Double_t Amplitude=2, Double_t landau_mean=1000,Double_t landau_sigma=1000,Double_t gauss_sigma=1000);
	void setLimits_Amplitude(Double_t startValue=0,Double_t endValue=2);
	void setLimits_LandauMP(Double_t startValue=0,Double_t endValue=1000);
	void setLimits_LandauSigma(Double_t startValue=0,Double_t endValue=1000);
	void setLimits_GaussSigma(Double_t startValue=0,Double_t endValue=1000);
	void setFitOptions(const char* options);
  void setFitRange(Double_t MinValue, Double_t MaxValue);
  void setLandauGauss_Separation(Double_t separation = 100);
	// the default values have to be implemented different in future 
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	// getters
	Double_t getAmplitude();
	Double_t getLandauMostProbable();
	Double_t getLandauSigma();
	Double_t getGaussSigma();
	Double_t getErrorOfAmplitude();
	Double_t getErrorOfLandauMP();
	Double_t getErrorOfLandauSigma();
	Double_t getErrorOfGaussSigma();
	Double_t getChiSqare();  /* since for now there are no errors supported for the individual measuring points one can not apply 
							    the function chi^2 = \sum (y_i - fit(x))^2/Sigma_i
							    therefor one has to use a different function to estimate the Goodness of the fit 
							    Chi^^2 = \sum (y_i - fit(x))^2/fit(x)
							    in cases where "y" is somehow related to an abundance of an occurrence of an incident, for example 
							    the firing of an detector, the value obtained by this method is is just a constant fraction of the "true chi square"
							    since "sigma^2 ~ fit(x)" */
							
	Double_t GetNDF();
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	// Workers
	// Fits the data contained in the TGraph object
	Int_t operator()(TGraph *fitData,Int_t firstBin=0,Int_t Lastbin=-1);
	// Fits the data contained in the TH1D object
  Int_t operator()(TH1D *fitData, Int_t firstBin = 0, Int_t Lastbin = -1);
	//Int_t fitLandau( TGraph* fitData );


	/* fast fits:
			they don't make copies of the Data set therefore they are not drawable */

	// Fits the data contained in the TGraph object
	Int_t FastFit(TGraph *fitData);
	
	// Fits the data contained in the TH1D object
	Int_t FastFit(TH1D *fitData);

	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	// display Methods

	// draws the Data set and the Fit function on top
	Int_t DrawfitFunction();

	// prints Out the results of the fits
	void printResults(std::ostream& out);
	void printResults();
	void saveFitToFile(const char* fileName);

	//////////////////////////////////////////////////////////////////////////
	TF1* getLandauGauss();
private:
	
#ifndef __CINT__  // there is no need to create a Dictionary for this since you cannot call it from outside the class anyway

	Size_t Size();
	void DeleteCopies();
	Double_t chiSqareTGraph();
	Double_t chiSqareTH1();
	Double_t R_square();
	Double_t  parLimitsLo[4];
	Double_t parLimitsHi[4];
	Double_t fiterrors[4];
	Double_t StartValues[4];
	Double_t parameters[4];
	Double_t chi;
	Double_t NDf;
 
	Bool_t DrawAble;  // this flag decides if one can draw the output or not. in FastFit() no copie of the datas is made therefore it is not Drawable.
	TGraph *g; // This TGraph will save the data points. 
	TF1 *ffit; // the data will be fitted with this function
	TF1* fit_Landau_gauss_;  //Landau Convoluted with a gauss
	//TF1 *fitLandau_;  //this are the Basic functions just for Comparison 
	std::string fitOptions_;
	TH1D* hData;
	int numOfFits;
#endif // __CINT__
	ClassDef(landgausFit, 1);
};

