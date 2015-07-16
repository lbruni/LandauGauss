#include "landgausFit.h"
#include "TF1.h"
#include <TH1D.h>
#include "TSpline.h"
#include <TMath.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <fstream>
#include "fitFunctionClass.hh"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#define lanMp   0
#define lanSig  1
#define Ampl    2
#define gausSig 3
#define parSize 4


using namespace std;
// first convolute than inegerating
Double_t LandauGaus(Double_t* x, Double_t* par);



// The Fit Method cannot access member functions. therefore this functions has to be global.
Double_t landau_function( Double_t* x, Double_t* para );
Double_t LandauGausInt( Double_t *x,Double_t *par );
void DrawCopyOfTGraph(TGraphErrors * graph);


TGraphErrors* makeCopieOfTH1D(TH1D* h1,Int_t firstBin=0,Int_t lastBin=-1);
TGraphErrors* makeCopieOfTGraph(TGraphErrors* h1, Int_t firstBin = 0, Int_t lastBin = -1);


landgausFit::landgausFit(void)
{
    m_c=nullptr;
    numOfFits=0;
    g=nullptr;
    hData=nullptr;
    DrawAble =false;
    setStartValues();
    
    ffit =nullptr;
    // preparing the fit function
    fit_Landau_gauss_ = new TF1("fitFun", &NewLandauGausInt, 0, 250, 4);
    fit_Landau_gauss_->SetParNames("landau_mean","landau_sigma","Amplitude","gaus_sigma");
    
    
    //preparing the Landau Function
    
    // setting everything to the Default Values
    setLimits_Amplitude();
    setLimits_GaussSigma();
    setLimits_LandauMP();
    setLimits_LandauSigma();
    setLandauGauss_Separation();
    fitOptions_="RB0Q";
}


landgausFit::~landgausFit(void)
{
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
    StartValues[lanMp]=landau_mean;
    StartValues[lanSig]=landau_sigma;
    StartValues[Ampl]=Amplitude;
    StartValues[gausSig]=gaus_sigma;
    
    
}

void landgausFit::setStartAmplitude( Double_t Amplitude/*=1*/ )
{
    StartValues[Ampl]=Amplitude;
}

void landgausFit::setStartLandauMean( Double_t landau_mean/*=200*/ )
{
    StartValues[lanMp]=landau_mean;
}

void landgausFit::setStartLandauSigma( Double_t landau_sigma/*=70*/ )
{
    StartValues[lanSig]=landau_sigma;
}

void landgausFit::setStartGaussSigma( Double_t gaus_sigma/*=10*/ )
{
    StartValues[gausSig]=gaus_sigma;
}

void landgausFit::setLowerLimits( Double_t Amplitude/*=0*/, Double_t landau_mean/*=0*/,Double_t landau_sigma/*=0*/,Double_t gaus_sigma/*=0*/ )
{
    parLimitsLo[lanMp]=landau_mean;
    parLimitsLo[lanSig]=landau_sigma;
    parLimitsLo[Ampl]=Amplitude;
    parLimitsLo[gausSig]=gaus_sigma;
}

void landgausFit::setUpperLimits( Double_t Amplitude/*=2*/, Double_t landau_mean/*=1000*/,Double_t landau_sigma/*=1000*/,Double_t gaus_sigma/*=1000*/ )
{
    parLimitsHi[lanMp]=landau_mean;
    parLimitsHi[lanSig]=landau_sigma;
    parLimitsHi[Ampl]=Amplitude;
    parLimitsHi[gausSig]=gaus_sigma;
}

void landgausFit::setLimits_Amplitude( Double_t startValue/*=0*/,Double_t endValue/*=2*/ )
{
    parLimitsLo[Ampl]=startValue;
    parLimitsHi[Ampl]=endValue;
}

void landgausFit::setLimits_LandauMP( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
    parLimitsLo[lanMp]=startValue;
    parLimitsHi[lanMp]=endValue;
}

void landgausFit::setLimits_LandauSigma( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
    parLimitsLo[lanSig]=startValue;
    parLimitsHi[lanSig]=endValue;
}

void landgausFit::setFitOptions( const char* options )
{
    fitOptions_=options;
}

void landgausFit::setLimits_GaussSigma( Double_t startValue/*=0*/,Double_t endValue/*=1000*/ )
{
    parLimitsLo[gausSig]=startValue;
    parLimitsHi[gausSig]=endValue;
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
    return parameters[Ampl];
}

Double_t landgausFit::getErrorOfAmplitude()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return fiterrors[Ampl];
}

Double_t landgausFit::getLandauSigma()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return parameters[lanSig];
}

Double_t landgausFit::getErrorOfGaussSigma()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return fiterrors[gausSig];
}

Double_t landgausFit::getGaussSigma()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return parameters[gausSig];
}

Double_t landgausFit::getErrorOfLandauSigma()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return fiterrors[lanSig];
}

Double_t landgausFit::getLandauMostProbable()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return parameters[lanMp];
}


Double_t landgausFit::getErrorOfLandauMP()
{
    if (numOfFits==0)
    {
        return ERROR_RETURN_VALUE;
    }
    return fiterrors[lanMp];
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

Int_t landgausFit::operator()(TGraphErrors *fitData, Int_t firstBin, Int_t Lastbin)
{
    
    // only one type of data container should be active therefore the program checks both container types
    DeleteCopies();
    g=makeCopieOfTGraph(fitData,firstBin,Lastbin);
    
    
    
    
    DrawAble=true;
    ++numOfFits;
    //////////////////////////////////////////////////////////////////////////
    // landgaus.C
    // Once again, here are the Landau * Gaussian parameters:
    //   par[lanMp]=Width (scale) parameter of Landau density
    //   par[lanSig]=Most Probable (MP, location) parameter of Landau density
    //   par[Ampl]=Total area (integral -inf to inf, normalization constant)
    //   par[gausSig]=Width (sigma) of convoluted Gaussian function
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
    
    ffit=fit_Landau_gauss_;
    ffit->SetParameters(StartValues);
    for (Int_t i=0; i<4; i++) {
        ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
    }
    
    TCanvas *canvas = new TCanvas("canvas", "", 500,500);
    g->Fit(ffit);//,fitOptions_.c_str());
    g->Draw("AP");
    
    canvas->SaveAs("plotinfitdovrebbeesserefittato.png");
    ffit->GetParameters(parameters);    // obtain fit parameters
    
    for (Int_t i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    chi = ffit->GetChisquare();  // obtain chi^2
    NDf = ffit->GetNDF();           // obtain ndf
    
    return SUCCESS_RETURN_VALUE;
}

Int_t landgausFit::operator()(TH1D *fitData, Int_t firstBin , Int_t Lastbin )
{
    
    DeleteCopies();
    g=makeCopieOfTH1D(fitData,firstBin,Lastbin);
    DrawAble=true;
    ++numOfFits;
    ffit=fit_Landau_gauss_;
    ffit->SetParameters(StartValues);
    
    
    for (Int_t i=0; i<4; i++) {
        ffit->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
    }
    
    g->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot
    
    ffit->GetParameters(parameters);    // obtain fit parameters
    for (Int_t i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    chi = ffit->GetChisquare();  // obtain chi^2
    NDf = ffit->GetNDF();           // obtain ndf
    
    //	return (ffit);              // return fit function
    
    
    return SUCCESS_RETURN_VALUE;
}
// Int_t landgausFit::fitLandau( TGraph* fitData )
// 	{
// 	DrawAble=true;
// 	++numOfFits;
// 	DeleteCopies();
//
// 	g=makeCopieOfTGraph(fitData);
//
//
// 	ffit=fitLandau_;
// 	fitLandau_->SetParameters(StartValues);
//
//
// 	for (Int_t i=0; i<3; i++) {
// 		fitLandau_->SetParLimits(i, parLimitsLo[i], parLimitsHi[i]);
// 	}
//
// 	fitData->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot
//
// 	fitLandau_->GetParameters(parameters);    // obtain fit parameters
// 	for (Int_t i=0; i<3; i++) {
// 		fiterrors[i] = fitLandau_->GetParError(i);     // obtain fit parameter errors
// 	}
// 	fiterrors[gausSig]=0; // the value should not be nothing
// 	parameters[gausSig]=0; // same here
//
// 	chi = fitLandau_->GetChisquare();  // obtain chi^2
// 	NDf = fitLandau_->GetNDF();           // obtain ndf
//
// 	//	return (ffit);              // return fit function
//
// 	return SUCCESS_RETURN_VALUE;
//
// }
Int_t landgausFit::FastFit( TGraphErrors *fitData )
{
    // only one type of data container should be active therefore the program checks both container types
    DeleteCopies();
    DrawAble =false;
    ++numOfFits;
    //////////////////////////////////////////////////////////////////////////
    // landgaus.C
    // Once again, here are the Landau * Gaussian parameters:
    //   par[lanMp]=Width (scale) parameter of Landau density
    //   par[lanSig]=Most Probable (MP, location) parameter of Landau density
    //   par[Ampl]=Total area (integral -inf to inf, normalization constant)
    //   par[gausSig]=Width (sigma) of convoluted Gaussian function
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
    
    g->Fit(ffit,fitOptions_.c_str());   // fit within specified range, use ParLimits, do not plot
    
    ffit->GetParameters(parameters);    // obtain fit parameters
    for (Int_t i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    chi = ffit->GetChisquare();  // obtain chi^2
    NDf = ffit->GetNDF();           // obtain ndf
    
    //	return (ffit);              // return fit function
    
    
    return SUCCESS_RETURN_VALUE;
}



Int_t landgausFit::DrawfitFunction(Int_t channel)
{
    TCanvas *c = new TCanvas();
    TFile *f = new TFile("output.root","UPDATE");
    
    if (DrawAble==false)
    {
        cout<<"Fast fit does not support Drawing the plot afterwards"<<endl;
        return ERROR_RETURN_VALUE;
    }
    if (g!=nullptr)
    {   //g->GetXaxis()->SetLimits(0,800);
        //g->GetYaxis()->SetLimits(0,1);
        g->Draw("AP");
    }
    if (hData!=nullptr)
    {
        hData->GetXaxis()->SetLimits(0,800);
        hData->GetYaxis()->SetLimits(0,1);
        hData->Draw("AB");
    }
    
    printResults();
    std::string outname = "Scurve_Channel_" + std::to_string(channel) + ".png";
    std::string graphname = "Scurve_Channel_" + std::to_string(channel);
    std::cout<<outname<<std::endl;
    g->GetXaxis()->SetTitle("Threshold [mV]");
    g->GetYaxis()->SetTitle("Efficiency");
    gStyle->SetOptStat();
    g->SetName(graphname.c_str());
    c->SaveAs(outname.c_str());
    
    TCanvas *c1 = new TCanvas();
    TF1 *f1;
    f1 = ffit;
    f1->SetParameter("Amplitude",-1);
    TGraph *g1 = (TGraph*)f1->DrawDerivative();
    std::string outname_derivative = "Derivative_Scurve_Channel_" + std::to_string(channel) + ".png";
    std::string graphname_derivative = "Derivative_Scurve_Channel_" + std::to_string(channel);
    g1->SetName(graphname_derivative.c_str());
    c1->SaveAs(outname_derivative.c_str());
    
    
    
    g->Write();
    g1->Write();
    f->Print();
    f->Close();
    
    
    return SUCCESS_RETURN_VALUE;
}

void landgausFit::printResults(ostream& out)
{
    
    
    out<<"Fit results: "<<endl;
    out<<"Landau Most Probable Value: "<<getLandauMostProbable()<<" +- "<< getErrorOfLandauMP() <<endl;
    out<<"Landau Sigma              : "<<getLandauSigma()<<" +- "<< getErrorOfLandauSigma()<<endl;
    out<<"Amplitude                 : "<<getAmplitude()<<" +- "<<getErrorOfAmplitude()<<endl;
    out<<"Gauss Sigma               : "<<getGaussSigma()<<" +- "<< getErrorOfGaussSigma()<<endl;
    out<<"Degrees of Freedom (NDF)  : "<<GetNDF()<<endl;
    out<<"Chi^2                     : "<<getChiSqare()<<endl;  // read the comments to the ChiSquare Methode
    out << "===============================================================================" << endl;
}



void landgausFit::printResults(){
    
    printResults(cout);
}

void landgausFit::saveFitToFile( const char* fileName )
{
    ofstream out(fileName);
    out << "Fit results: " << endl;
    out<<"Landau Most Probable Value: "<<getLandauMostProbable()<<" +- "<< getErrorOfLandauMP() << endl;
    out<<"Landau Sigma              : "<<getLandauSigma()<<" +- "<<getErrorOfLandauSigma()<<endl;
    out<<"Amplitude                 : "<<getAmplitude() <<" +- "<<getErrorOfAmplitude()<<endl;
    out<<"Gauss Sigma               : "<< getGaussSigma()<<" +- "<<getErrorOfGaussSigma()<<endl;
    out<<"data: "<<endl;
    out << "x ;  data_y ; Fit_y" << endl;
    
    
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

TF1* landgausFit::getLandauGauss()
{
    return fit_Landau_gauss_;
}

void landgausFit::setFitRange(Double_t MinValue, Double_t MaxValue)
{
    fit_Landau_gauss_->SetRange(MinValue, MaxValue);
}

void landgausFit::setLandauGauss_Separation(Double_t separation /*= 100*/)
{
    SetNewLandauGauss_Setparation(separation);
}





TGraphErrors* makeCopieOfTH1D(TH1D* h1, Int_t firstBin/*=0*/, Int_t lastBin/*=-1*/)
{
    // somehow i had some problems with just 
    // cloning the data therefore i make the copy by hand
    
    TGraphErrors *g1 = new TGraphErrors();
    if (lastBin==-1)
    {
        lastBin = h1->GetNbinsX();
    }
    Int_t counter = 0;
    for (size_t i = firstBin; i < lastBin; ++i)
    {
        g1->SetPoint(counter++, h1->GetBinCenter(i), h1->GetBinContent(i));
    }
    
    g1->SetEditable(false);
    
    return g1;
}


TGraphErrors* makeCopieOfTGraph(TGraphErrors* graph_in, Int_t firstBin /*= 0*/, Int_t lastBin /*= -1*/)
{
    // somehow i had some problems with just 
    // cloning the data therefore i make the copy by hand
    TGraphErrors *g1 = new TGraphErrors();
    if (lastBin == -1)
    {
        lastBin = graph_in->GetN();
    }
    
    Int_t counter = 0;
    for (size_t i = firstBin; i <lastBin; ++i)
    {
        Double_t x, y, ex, ey;
        graph_in->GetPoint(i, x, y);
        ex= graph_in->GetErrorX(i);//lb
        ey= graph_in->GetErrorY(i);//lb
        g1->SetPoint(i, x, y);
        g1->SetPointError(i, ex, ey);//lb
    }
    
    g1->SetEditable(false);
    DrawCopyOfTGraph(g1);
    return g1;
}

void DrawCopyOfTGraph(TGraphErrors * graph)
{
    TCanvas *c = new TCanvas("c","",500, 500);
    graph->Draw("AP");
    c->SaveAs("provagraphinlib.png");
    
}