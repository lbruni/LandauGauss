#include <windows.h>
#include <tchar.h> 
#include <stdio.h>
#include <strsafe.h>
#include <vector>
#include "TFile.h"
#include <TH2.h>
#include <TH1.h>
#include <TObjString.h>
#include <iostream>
#include <fstream>


#include <TGraph.h>

#define SINGLEEVENT

#include "landgausFit.h"

#define ReadOutChannels 2560
#define bufferSize 20000
std::vector<WIN32_FIND_DATA> dir(const TCHAR* param);
TH1D * get_strip_effi(TH2* inputHist);
char * WChar2Char( wchar_t *orig);
void DisplayErrorBox(LPTSTR lpszFunction);

#pragma comment(lib, "User32.lib")
using namespace std;



TH2D loopOverFolder( Int_t &size)
{

	

	auto files=dir(TEXT("*th*.root"));

	TH2D effi("strip_Effi","Effi VS strip And threshold",ReadOutChannels,0,ReadOutChannels,files.size(),0,files.size());
	//Threshold=new Double_t[files.size()];
	//vector<Double_t> threshold;
	TH1D* h1;
	Int_t j=1;
	for (auto e:files)
	{
		_tprintf(TEXT("Current File %s\n"),e.cFileName);


		TFile inputFile((e.cFileName));
		TH2D *h2=(TH2D*)inputFile.Get("TH2;1");
		h1=get_strip_effi(h2);
		TObjString *s=(TObjString*)inputFile.Get("Threshold;1");
		//threshold.push_back(stod(s->GetName()));
		//Threshold[j-1]=stod(s->GetName());
		//cout<<Threshold[j-1]<<endl;
		effi.SetBinContent(0,j,stod(s->GetName()));
		for (int k=1;k<ReadOutChannels;++k)
		{
			effi.SetBinContent(k,j,h1->GetBinContent(k));
		}
		++j;

		//threshold.push_back(stod(s->GetName()));

	}
	size=files.size();
	return std::move(effi);
}


TGraph getStrip(TH2D *effi,Int_t stripNumber){

	//TGraph g2;
    
	TH1D* h11=effi->ProjectionY("asddasd",stripNumber,stripNumber);
	TGraph g1;
	//Double_t numberOfEvents=0;
	for (size_t i=0;i<h11->GetNbinsX()-1;++i)
	{
		g1.SetPoint(i,effi->GetBinContent(0,i+1),h11->GetBinContent(i+1)); // in this project the underflow bin of the 2d histogramm is used as X axis
	//	cout<<effi->GetBinContent(0,i+1)<<"        "<<h11->GetBinContent(i+1)<<endl;
	//	numberOfEvents+=h11->GetBinContent(i+1);


	}
	g1.SetEditable(false);

	return std::move(g1);
}

TGraph Analyse(TH2D* effi, Int_t size){
	TGraph g2;

	Double_t *Threshold=new Double_t[size];
	landgausFit fit1;
	for (int i = 0; i < size; i++)
	{
		Threshold[i]=effi->GetBinContent(0,i+1);
	//	cout<<Threshold[i]<<endl;
	}
	

	for(int j=200;j<300;++j){


	TH1D* h11=effi->ProjectionY("asddasd",j,j);
	TGraph g1;
	Double_t numberOfEvents=0;
	for (size_t i=0;i<size;++i)
	{
		g1.SetPoint(i,Threshold[i],h11->GetBinContent(i+1));
		numberOfEvents+=h11->GetBinContent(i+1);

	}
	g1.SetEditable(false);

	if (numberOfEvents>0.01)
	{

		fit1(&g1);




		fit1.printResults();
		g2.SetPoint(j,j,fit1.getAmplitude());
	}else

	{
		g2.SetPoint(j,j,0);
	}
		
}

//g2.Draw("aL*");


	return std::move(g2);
	
}

std::vector<WIN32_FIND_DATA> dir(const TCHAR* param){
	std::vector<WIN32_FIND_DATA> returnValue;

	WIN32_FIND_DATA Filebuffer;


	HANDLE hFind = INVALID_HANDLE_VALUE;
	hFind = FindFirstFile(param, &Filebuffer);
	if (INVALID_HANDLE_VALUE == hFind) 
	{
		_tprintf(TEXT("error"));
	} else 
	{


		do
		{
			returnValue.push_back(std::move(Filebuffer));
		}
		while (FindNextFile(hFind, &Filebuffer) != 0);
	}
	return std::move(returnValue);

}


char * WChar2Char( wchar_t *orig){
	// Convert to a char*
	size_t origsize = wcslen(orig) + 1;
	const size_t newsize = 100;
	size_t convertedChars = 0;
	char *nstring=new char[newsize];
	wcstombs_s(&convertedChars, nstring, origsize, orig, _TRUNCATE);
	//strcat_s(nstring, " (char *)");
	cout << nstring << endl;
	return nstring;

}

TH1D * get_strip_effi(TH2* inputHist){

	TH1D *returnValue=inputHist->ProjectionX("hitMap");
	Int_t NumOfRows=inputHist->GetNbinsY();
	for (Int_t i=1;i<returnValue->GetNbinsX();++i)
	{
		returnValue->SetBinContent(i,returnValue->GetBinContent(i)/NumOfRows);
	}
	return returnValue;
}
Double_t cutTheString(std::string &inputstr){
	Double_t returnValue=-1;
	if (inputstr.size()>0)
	{
		int index=inputstr.find_first_not_of("0123456789+-. ");
		if (index>0)
		{
			returnValue=std::stod(inputstr.substr(0,index));
			inputstr=inputstr.substr(index+1);

		}

	}

	return returnValue;

}

Int_t getSizeOfFile(const char * fileName){
	Int_t returnValue=0;

	ifstream inFile(fileName);
	char buffer[bufferSize];
	inFile.getline(buffer,bufferSize); // first line contains header information.
	while(!(inFile.eof())){
		inFile.getline(buffer,bufferSize);

		returnValue++;
		if (returnValue>10000)
		{returnValue-1;
			break;
		}
	}

	std::string s(buffer);
	if (s.size()<2&&returnValue>0)
	{
		returnValue--;
	}
	return returnValue;
}

TH2D* load_ascii2(const char* fileName,int x_bins,Double_t &Threshold){

	Int_t sizeOfFile=getSizeOfFile(fileName);
	ifstream inFile(fileName);
	std::string buffer;
	getline(inFile,buffer);
#ifdef _DEBUG
	std::cout<<buffer<<std::endl;
#endif // _DEBUG


	char HistogramTitel[20];
	Threshold=std::stod(buffer);
	sprintf(HistogramTitel,"Threshold %s mV",buffer.c_str());
	TH2D *eventbuffer=new TH2D("TH2",HistogramTitel,x_bins,0,x_bins,sizeOfFile,0,sizeOfFile);
	Int_t x=0,y=0;
	Double_t value1;
	while (!inFile.eof())
	{
		try
		{


			getline(inFile,buffer,';');
			value1=std::stod(buffer);
			//std::cout<<value1;
			eventbuffer->Fill(x++,y,value1);
			if (x>=x_bins)
			{
				x=0;
				++y;
				//std::cout<<std::endl;
			}
		}
		catch (...){
#ifdef  _DEBUG
			std::cout<<"error"<<std::endl;
#endif //  _DEBUG



		}
	}
	return eventbuffer;

}

TH2D load_ascii3(const char* filename){

	Int_t sizeOfFile=getSizeOfFile(filename);
//	cout<<sizeOfFile<<endl;
	ifstream inFile(filename);
	//char buffer[bufferSize];
	std::string header;
	getline(inFile,header);  // Gives the header


	std::string FirstLine;
	getline(inFile,FirstLine);  // from the first line it estimated the number of columns



	size_t x_bins = std::count(FirstLine.begin(), FirstLine.end(), ';');
//	cout<<x_bins<<endl;
	TH2D eventbuffer("TH2",header.c_str(),x_bins,0,x_bins,sizeOfFile,0,sizeOfFile);
	std::string buffer;


	int i=0;
	Double_t value1=0;
	Int_t x=0,y=0;
	while (FirstLine.size()>0&&value1>=0)
	{
		value1=cutTheString(FirstLine);
		eventbuffer.SetBinContent(i++,y,value1);

	}



	value1=0;

	++y;



	while (!inFile.eof())
	{
		try
		{



			getline(inFile,buffer,';');
			value1=std::stod(buffer);
			//std::cout<<value1;
			eventbuffer.SetBinContent(x++,y,value1);
			if (x>=x_bins)
			{
				x=0;
				++y;
				//std::cout<<std::endl;
			}
		}
		catch (...){
#ifdef  _DEBUG
			std::cout<<"error"<<std::endl;
#endif //  _DEBUG



		}


	}
	return std::move(eventbuffer);
}

void save2ascii(const char* name, const TH2D& inputHist ){

	ofstream out(name);
for (Int_t i=0; i<inputHist.GetNbinsY();++i)
{

	for (Int_t j=0; j<inputHist.GetNbinsX();++j)
	{
		out<<inputHist.GetBinContent(j,i)<<" ; ";
	}
	out<<endl;
}
out.close();
}