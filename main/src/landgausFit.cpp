#include "landgausFit.h"
#include "TF1.h"
#include <TH1D.h>
#include "TSpline.h"
#include <TMath.h>
//#include <TAxis.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

// The Fit Method cannot access member functions. therefore this functions has to be global.
	Double_t landau_function( Double_t* x, Double_t* para );
	Double_t LandauGaus( Double_t *x,Double_t *par );
	Double_t landau_function_int( Double_t* x, Double_t* para );
	
	// creates a spline Function of the the landau Function. this function will be called Automatically by the creation of the object. The spline will only be created independing on how many instances of this class you have
	void makeSplineFunction(const Double_t StartOfTheInterval=-100,const Double_t endOfTheInterval=500,const size_t steps=300);
	TSpline3* g_spline=nullptr;
	Double_t gSplineMin=-10,gSplineMax=40;
//

	TGraph* makeCopieOfTH1D(TH1D* h1);

landgausFit::landgausFit(void)
{
	numOfFits=0;
	g=nullptr;
	hData=nullptr;
	DrawAble =false;
// 	parameters=new Double_t[4];
// 	StartValues=new Double_t[4];
// 	fiterrors=new Double_t[4];
	setStartValues();


	makeSplineFunction();//-100,500,300);

	ffit =nullptr;
	// preparing the fit function
	fit_Landau_gauss_ = new TF1("fitFun",&LandauGaus,0,1000,4);
	fit_Landau_gauss_->SetParNames("landau_mean","landau_sigma","Amplitude","gaus_sigma");

	
	//preparing the Landau Function
	fitLandau_=new TF1("landau_int",&landau_function_int,0,1000,3);
	fitLandau_->SetParNames("landau_mean","landau_sigma","Amplitude");
	
	// setting everything to the Default Values
	setLimits_Amplitude();
	setLimits_GaussSigma();
	setLimits_LandauMP();
	setLimits_LandauSigma();
	fitOptions_="RB0Q";
}


landgausFit::~landgausFit(void)
{
	
	if (g_spline!=nullptr)
	{
		delete g_spline;
		g_spline=nullptr;

	}

	if (g!=nullptr)
	{
		delete g;
		g=nullptr;
	}
	if (hData!=nullptr)
	{
		delete hData;
		hData =nullptr;
	}
	if (fit_Landau_gauss_!=nullptr)
	{
		delete fit_Landau_gauss_;
		fit_Landau_gauss_=nullptr;
	}
	if (fitLandau_!=nullptr)
	{
		delete fitLandau_;
		fitLandau_=nullptr;
	}
	
}
void landgausFit::DeleteCopies()
{
	// only one type of data container should be active therefore the program checks both container types 
	if (g!=nullptr)
	{
		delete g;
		g=nullptr;
	}
	if (hData!=nullptr)
	{
		delete hData;
		hData=nullptr;
	}
}

//////////////////////////////////////////////////////////////////////////
// Setter Functions
//////////////////////////////////////////////////////////////////////////

void landgausFit::setStartValues( Double_t Amplitude/*=1*/, Double_t landau_mean/*=200*/,Double_t landau_sigma/*=70*/,Double_t gaus_sigma/*=10*/ )
{
	StartValues[0]=landau_mean;
	StartValues[1]=landau_sigma;
	StartValues[2]=Amplitude;
	StartValues[3]=gaus_sigma;


}

void landgausFit::setStartAmplitude( Double_t Amplitude/*=1*/ )
{
	StartValues[2]=Amplitude;
}

void landgausFit::setStartLandauMean( Double_t landau_mean/*=200*/ )
{
		StartValues[0]=landau_mean;
}

void landgausFit::setStartLandauSigma( Double_t landau_sigma/*=70*/ )
{
		StartValues[1]=landau_sigma;
}

void landgausFit::setStartGaussSigma( Double_t gaus_sigma/*=10*/ )
{
	StartValues[3]=gaus_sigma;
}

void landgausFit::setLowerLimits( Double_t Amplitude/*=0*/, Double_t landau_mean/*=0*/,Double_t landau_sigma/*=0*/,Double_t gaus_sigma/*=0*/ )
{
	parLimitsLo[0]=landau_mean;
	parLimitsLo[1]=landau_sigma;
	parLimitsLo[2]=Amplitude;
	parLimitsLo[3]=gaus_sigma;
}

void landgausFit::setUpperLimits( Double_t Amplitude/*=2*/, Double_t landau_mean/*=1000*/,Double_t landau_sigma/*=1000*/,Double_t gaus_sigma/*=1000*/ )
{
	parLimitsHi[0]=landau_mean;
	parLimitsHi[1]=landau_sigma;
	parLimitsHi[2]=Amplitude;
	parLimitsHi[3]=gaus_sigma;
}

void landgausFit::setLimits_Amplitude( Double_t startValue/*=0*/,Double_t endValue/*=2*/ )
{
	parLimitsLo[2]=startValue;
	parLimitsHi[2]=endValue;
}

void landgausFit::setLimits_LandauMP( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
	parLimitsLo[0]=startValue;
	parLimitsHi[0]=endValue;
}

void landgausFit::setLimits_LandauSigma( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
	parLimitsLo[1]=startValue;
	parLimitsHi[1]=endValue;
}

void landgausFit::setFitOptions( const char* options )
{
	fitOptions_=options;
}

void landgausFit::setLimits_GaussSigma( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
	parLimitsLo[3]=startValue;
	parLimitsHi[3]=endValue;
}
//////////////////////////////////////////////////////////////////////////
// end Setter Functions
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
// Getter Functions 
//////////////////////////////////////////////////////////////////////////
Double_t landgausFit::getAmplitude()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return parameters[2];
}

Double_t landgausFit::getErrorOfAmplitude()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return fiterrors[2];
}

Double_t landgausFit::getLandauSigma()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return parameters[1];
}

Double_t landgausFit::getErrorOfGaussSigma()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return fiterrors[3];
}

Double_t landgausFit::getGaussSigma()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return parameters[3];
}

Double_t landgausFit::getErrorOfLandauSigma()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return fiterrors[1];
}

Double_t landgausFit::getLandauMostProbable()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return parameters[0];
}


Double_t landgausFit::getErrorOfLandauMP()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return fiterrors[0];
}

Double_t landgausFit::GetNDF()
{
	if (numOfFits==0)
	{
		return ERROR_RETURN_VALUE;
	}
	return NDf;
}

//////////////////////////////////////////////////////////////////////////
// end Getter Functions
//////////////////////////////////////////////////////////////////////////




Int_t landgausFit::operator()( TGraph *fitData )
{

	// only one type of data container should be active therefore the program checks both container types 
	DeleteCopies();
	g=dynamic_cast<TGraph*>(fitData->Clone());
	DrawAble=true;
	++numOfFits;
	//////////////////////////////////////////////////////////////////////////
	// landgaus.C
	// Once again, here are the Landau * Gaussian parameters:
	//   par[0]=Width (scale) parameter of Landau density
	//   par[1]=Most Probable (MP, location) parameter of Landau density
	//   par[2]=Total area (integral -inf to inf, normalization constant)
	//   par[3]=Width (sigma) of convoluted Gaussian function
	//
	// Variables for langaufit call:
	//   his             histogram to fit
	//   fitrange[2]     lo and hi boundaries of fit range
	//   startvalues[4]  reasonable start values for the fit
	//   parlimitslo[4]  lower parameter limits
	//   parlimitshi[4]  upper parameter limits
	//   fitparams[4]    returns the final fit parameters
	//   fiterrors[4]    returns the final fit errors
	//   ChiSqr          returns the chi square
	//   NDF             returns ndf

	
	

	//sprintf(FunName,"Fitfcn_%s",his->GetName());

	ffit=fit_Landau_gauss_;

	ffit->SetParameters(StartValues);

	//parLimitsLo={0,0,0,0};
	 //Double_t parLimitsHi[4]={1000,1000,2,1000};
	for (Int_t i=0; i<4; i++) {
		ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
	}

	//fitData->Fit(ffit,"RB0QWW");   // fit within specified range, use ParLimits, do not plot
	fitData->Fit(ffit,fitOptions_.c_str());   
	ffit->GetParameters(parameters);    // obtain fit parameters
	for (Int_t i=0; i<4; i++) {
		fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
	}
	 chi = ffit->GetChisquare();  // obtain chi^2
	 NDf = ffit->GetNDF();           // obtain ndf
	 
//	return (ffit);              // return fit function

return SUCCESS_RETURN_VALUE;
}

Int_t landgausFit::operator()( TH1D *fitData )
{
	
DeleteCopies();
g=makeCopieOfTH1D(fitData);
	DrawAble=true;
	++numOfFits;
	ffit=fit_Landau_gauss_;
	ffit->SetParameters(StartValues);


	for (Int_t i=0; i<4; i++) {
		ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
	}

	fitData->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot

	ffit->GetParameters(parameters);    // obtain fit parameters
	for (Int_t i=0; i<4; i++) {
		fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
	}
	chi = ffit->GetChisquare();  // obtain chi^2
	NDf = ffit->GetNDF();           // obtain ndf

	//	return (ffit);              // return fit function

	
	return SUCCESS_RETURN_VALUE;
}
Int_t landgausFit::fitLandau( TGraph* fitData )
	{
	DrawAble=true;
	++numOfFits;
	DeleteCopies();

	g=dynamic_cast<TGraph*>(fitData->Clone());

	
	ffit=fitLandau_;
	fitLandau_->SetParameters(StartValues);


	for (Int_t i=0; i<3; i++) {
		fitLandau_->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
	}

	fitData->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot

	fitLandau_->GetParameters(parameters);    // obtain fit parameters
	for (Int_t i=0; i<3; i++) {
		fiterrors[i] = fitLandau_->GetParError(i);     // obtain fit parameter errors
	}
	fiterrors[3]=0; // the value should not be nothing
	parameters[3]=0; // same here

	chi = fitLandau_->GetChisquare();  // obtain chi^2
	NDf = fitLandau_->GetNDF();           // obtain ndf

	//	return (ffit);              // return fit function

	return SUCCESS_RETURN_VALUE;

}
Int_t landgausFit::FastFit( TGraph *fitData )
	{
	// only one type of data container should be active therefore the program checks both container types 
	DeleteCopies();
	DrawAble =false;
	++numOfFits;
	//////////////////////////////////////////////////////////////////////////
	// landgaus.C
	// Once again, here are the Landau * Gaussian parameters:
	//   par[0]=Width (scale) parameter of Landau density
	//   par[1]=Most Probable (MP, location) parameter of Landau density
	//   par[2]=Total area (integral -inf to inf, normalization constant)
	//   par[3]=Width (sigma) of convoluted Gaussian function
	//
	// Variables for langaufit call:
	//   his             histogram to fit
	//   fitrange[2]     lo and hi boundaries of fit range
	//   startvalues[4]  reasonable start values for the fit
	//   parlimitslo[4]  lower parameter limits
	//   parlimitshi[4]  upper parameter limits
	//   fitparams[4]    returns the final fit parameters
	//   fiterrors[4]    returns the final fit errors
	//   ChiSqr          returns the chi square
	//   NDF             returns ndf




	//sprintf(FunName,"Fitfcn_%s",his->GetName());

	ffit=fit_Landau_gauss_;

	ffit->SetParameters(StartValues);

	//parLimitsLo={0,0,0,0};
	//Double_t parLimitsHi[4]={1000,1000,2,1000};
	for (Int_t i=0; i<4; i++) {
		ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
	}

	fitData->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot

	ffit->GetParameters(parameters);    // obtain fit parameters
	for (Int_t i=0; i<4; i++) {
		fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
	}
	chi = ffit->GetChisquare();  // obtain chi^2
	NDf = ffit->GetNDF();           // obtain ndf

	

	return SUCCESS_RETURN_VALUE;
	}

Int_t landgausFit::FastFit( TH1D *fitData )
{

	DeleteCopies();
	g=makeCopieOfTH1D(fitData);
	DrawAble=false;
	++numOfFits;
	ffit=fit_Landau_gauss_;
	ffit->SetParameters(StartValues);


	for (Int_t i=0; i<4; i++) {
		ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
	}

	fitData->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot

	ffit->GetParameters(parameters);    // obtain fit parameters
	for (Int_t i=0; i<4; i++) {
		fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
	}
	chi = ffit->GetChisquare();  // obtain chi^2
	NDf = ffit->GetNDF();           // obtain ndf

	//	return (ffit);              // return fit function


	return SUCCESS_RETURN_VALUE;
}



Int_t landgausFit::DrawfitFunction()
{
	if (DrawAble==false)
	{
		cout<<"Fast fit does not support Drawing the plot afterwards"<<endl;
		return ERROR_RETURN_VALUE;
	}
if (g!=nullptr)
{

	g->GetXaxis()->SetLimits(0,800);
	g->GetYaxis()->SetLimits(0,1);
	g->Draw("AP*");
}
if (hData!=nullptr)
{
	hData->GetXaxis()->SetLimits(0,800);
	hData->GetYaxis()->SetLimits(0,1);
	hData->Draw("AB");
}
	ffit->Draw("same");

	return SUCCESS_RETURN_VALUE;
}

void landgausFit::printResults(ostream& out)
{

	out<<"==============================================================================="<<endl;
	out<<"Fit results: "<<endl;
	out<<"Landau Most Probable Value: "<<getLandauMostProbable()<<" +- "<< getErrorOfLandauMP() <<endl;
	out<<"Landau Sigma              : "<<getLandauSigma()<<" +- "<< getErrorOfLandauSigma()<<endl;
	out<<"Amplitude                 : "<<getAmplitude()<<" +- "<<getErrorOfAmplitude()<<endl;
	out<<"Gauss Sigma               : "<<getGaussSigma()<<" +- "<< getErrorOfGaussSigma()<<endl;
	out<<"Degrees of Freedom (NDF)  : "<<GetNDF()<<endl;
	out<<"Chi^2                     : "<<getChiSqare()<<endl;  // read the comments to the ChiSquare Methode
}



void landgausFit::printResults(){

	printResults(cout);
}

void landgausFit::saveFitToFile( const char* fileName )
{
	ofstream out(fileName);
	out<<"Fit results: "<<endl;
	out<<"Landau Most Probable Value: "<<endl;
	out<<getLandauMostProbable()<<endl;
	out<<" +- "<<endl;
	out<<getErrorOfLandauMP() <<endl;

	out<<"Landau Sigma              : "<<endl;
	out<<getLandauSigma()<<endl;
	out<<" +- "<<endl;
	out<<getErrorOfLandauSigma()<<endl;
	out<<"Amplitude                 : "<<endl;
	out<<getAmplitude()<<endl;
	out<<" +- "<<endl;
	out<<getErrorOfAmplitude()<<endl;
	out<<"Gauss Sigma               : "<<endl;
	out<<getGaussSigma()<<endl;
	out<<" +- "<<endl;
	out<<getErrorOfGaussSigma()<<endl;
	out<<"data: "<<endl;

	Double_t x=-1;
	Double_t y=-1;
	Double_t y_fit=-1;

	for (Size_t i=0;i<Size();++i)
	{
	if (g!=nullptr)
	{
		g->GetPoint(i,x,y);
	}else if (hData!=nullptr)
	{
		x=hData->GetBinCenter(i);
		y=hData->GetBinContent(i);
	}


	if (ffit!=nullptr)
	{
		y_fit=ffit->Eval(x);
	}
		
		out<<x<<" ; "<< y<<"  ;  "<< y_fit <<endl;
	}

}

Double_t landgausFit::getChiSqare()
{
	Double_t returnValue=ERROR_RETURN_VALUE;
	if (g!=nullptr)
	{
		returnValue=chiSqareTGraph();

	}
	if (hData!=nullptr)
	{
		returnValue = chiSqareTH1();
	}
	return returnValue;
}

Double_t landgausFit::chiSqareTGraph()
{
	Double_t returnValue=0; //chiSquare

	for (Int_t i=0;i<g->GetN();++i)
	{
		Double_t x,y;
		g->GetPoint(i,x,y);
	//	cout<<x<<"  "<< y<<endl;
		if (ffit->Eval(x)>0)
		{
			returnValue+=(y-ffit->Eval(x))*(y-ffit->Eval(x))/ffit->Eval(x);
		}else{
			returnValue+=(y-ffit->Eval(x))*(y-ffit->Eval(x));
		}
		
	}
	return returnValue;
}

Double_t landgausFit::chiSqareTH1()
{
	Double_t returnValue=0; //chiSquare
	for (Int_t i=1;i<hData->GetNbinsX();++i)
	{
		Double_t x,y;
		x=hData->GetBinCenter(i);
		y=hData->GetBinContent(i);
		if (ffit->Eval(x)>0)
		{
		returnValue+=(y-ffit->Eval(x))*(y-ffit->Eval(x))/ffit->Eval(x);
		}else
		{
			returnValue+=(y-ffit->Eval(x))*(y-ffit->Eval(x));
		}

	}
	return returnValue;
}

Size_t landgausFit::Size()
{
	Double_t returnValue=ERROR_RETURN_VALUE;
	if (g!=nullptr)
	{
		returnValue=g->GetN();

	}
	if (hData!=nullptr)
	{
		returnValue = hData->GetNbinsX();
	}
	return returnValue;
}

void landgausFit::setFitRange(Double_t MinValue, Double_t MaxValue)
{
  fitLandau_->SetRange(MinValue, MaxValue);
}

void makeSplineFunction(const Double_t StartOfTheInterval/*=-10*/,const Double_t endOfTheInterval/*=40*/,const size_t steps/*=200*/)
{
	if (g_spline==nullptr)
	{
	gSplineMin=StartOfTheInterval;
	gSplineMax=endOfTheInterval;
	Double_t stepSize=(endOfTheInterval-StartOfTheInterval)/(steps-1);
	//////////////////////////////////////////////////////////////////////////
	// Susannes Code
	//
	const Double_t epsilon=1e-6;

	Double_t* x=new Double_t[steps];
	Double_t* y=new Double_t[steps];


	y[0]=1; // define integral as 0 for x=-n_sigmas*width
	x[0]=StartOfTheInterval;
	TF1 f1("lan","TMath::Landau(x,0,1,0)",StartOfTheInterval,endOfTheInterval);
	//f1.SetParameter(0,sigma);
	for(Int_t i=1;i<steps;++i) {
		x[i]=x[i-1]+stepSize;
		Double_t x1=x[i-1];
		Double_t x2=x[i];
		y[i]=y[i-1]-f1.Integral(x1,x2,(const Double_t*)0,epsilon);
	}

	g_spline=new TSpline3("Ilandau",x,y,steps,"b2e2",0,0);
	delete [] x;
	delete [] y;
	//////////////////////////////////////////////////////////////////////////
	}
}

Double_t landau_function( Double_t* x, Double_t* para )
{
	Double_t returnValue;
	//////////////////////////////////////////////////////////////////////////
	// Susannes Code
	//para 0=mpv
	//para 1=sigma
	Double_t dsig=((x[0]-para[0])/para[1]); //calculate distance from mpv in sigmas

	// a spline function only gives reasonable values inside the range it is defined
	// therefore one has to make sure that the values outside this range are reasonable

	if (dsig > gSplineMin && dsig < gSplineMax)
	{
		returnValue=g_spline->Eval(dsig);
	}
	else if (dsig <= gSplineMin)
	{
		returnValue=1;
	}else if (dsig>=gSplineMax)
	{
		returnValue=0;
	}

	return returnValue;
	//////////////////////////////////////////////////////////////////////////
}

Double_t LandauGaus( Double_t *x,Double_t *par )
{
	//////////////////////////////////////////////////////////////////////////
	// Susannes Code
	  TString funcstr="0";
	  Double_t tstep=par[3]/20;
      Int_t nsig=3;
  Double_t normgaus=0;
  
  //Zur Normierung: Schleife ueber nsig*sigma von Gauss (z.B. 2sigma). Um Summe des Gauss (ueber +-2sigma) um 0 zu berechnen. 
/*  
  for(Double_t t=-nsig*par[3];t<=nsig*par[3];t+=tstep) {
    normgaus+=TMath::Gaus(t,0,par[3]);
  }
*/  
  Double_t faltung=0;

  for(Double_t t=-nsig*par[3];t<=nsig*par[3];t+=tstep) {
    Double_t gaus=TMath::Gaus(-t,0,par[3],kTRUE);
    normgaus+=gaus;    
    Double_t xlandau=x[0]+t;
    faltung+=gaus*landau_function(&xlandau,par);
  }
 // cout<<" normgaus = "<<normgaus<<"   1/tstep = "<<1/tstep<<endl;
  faltung*=tstep;
  //faltung/=normgaus;
  faltung=faltung*par[2];



/*
  
  //Faltung durch Addition:Schleife ueber nsig*sigma von Gauss (z.B. 2sigma).
  // Um Summe des normierten Gauss (ueber +-2sigma) um 0 und des integrierten Landau an der Stelle x+t zu berechnen.
  for(Double_t t=-nsig*par[3];t<=nsig*par[3];t+=tstep) {
    funcstr+=(func_gspline(x+t,par[0],par[1]),TMath::Gaus(t,0,par[3])/normgaus,t);
  }*/

//////////////////////////////////////////////////////////////////////////
  
  return faltung;
}

Double_t landau_function_int( Double_t* x, Double_t* para ){
	//para 0=mpv
	//para 1=sigma
	//para 2=Amplitude
	return landau_function(x,para)*para[2];

}

TGraph* makeCopieOfTH1D(TH1D* h1){ // somehow i had some problems with just 
								   // cloning the data therefore i make the copy by hand

	TGraph *g1=new TGraph();
	for (size_t i=1;i<h1->GetNbinsX();++i)
	{
		g1->SetPoint(i-1,h1->GetBinCenter(i),h1->GetBinContent(i)); 
		Double_t x,y;
		g1->GetPoint(i-1,x,y);
		//cout<<h1->GetBinCenter(i)<<"   ;   " <<x<< "   ;   "<<h1->GetBinContent(i)<<"  ;  "<<y<<endl;
	}

	g1->SetEditable(false);

return g1;
}