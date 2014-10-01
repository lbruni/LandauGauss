#include "fitFunctionClass.hh"
#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <math.h>
#include <assert.h>
#include "TSpline.h"



using namespace std;





#define printVector(vec) print2file(#vec ## ".txt",vec)
template<typename T>
void print2file(const char* fileName, T t){
  std::ofstream out(fileName);

  for (auto&e : t){
    out << e << std::endl;
  }

}

axis makeLine(value_t start, value_t delta, value_t stop){
  

  auto d = stop - start;
  size_t vector_size = static_cast<Int_t>(round((stop - start) / delta)) + 1;

  axis ret(vector_size);

  value_t old = start - delta;
  for (auto&e:ret)
  {
    e = old + delta;
    old = e;

  }
  return ret;

}


template<typename Function,typename Container, typename... ARGGs>
Container callFunction(Function f, Container x_container, ARGGs... args){

  
  Container y_container(x_container.size());
  for (size_t i = 0; i < x_container.size(); ++i)
  {
    y_container[i] = f(x_container[i], args...,1);
  }

  return y_container;
}

template <typename Container>
Container convolut(Container A_Container, Container B_container){

  Container conv(A_Container.size() + B_container.size() + 1);

  for (int k = 0; k < conv.size();++k)
  {
    Container::value_type c = 0;

    for (int j = 0; j < A_Container.size();++j)
    {
      if (k + 1 - j>0 && k + 1 - j< B_container.size())
      {
        c += A_Container[j] * B_container[k + 1 - j];
      }

    }

    conv[k] = c;
  }

  return conv;
}


template <typename Container_t,typename VALUE_T>
Container_t cumsum(Container_t inContainer, VALUE_T startValue=0, VALUE_T Factor=1){
  Container_t ret(inContainer.size());

    ret[0] = startValue;
  for (size_t i = 1; i < ret.size(); ++i)
   {
      ret[i] = ret[i - 1] +Factor * inContainer[i - 1];
   }

    return ret;
}

Double_t NewLandauGausInt(Double_t *x, Double_t *par)
{
  static fitFunctionClass fitfun;
  fitfun.setNewParameter(par);
  return fitfun(x[0]);
}

void fitFunctionClass::newFunction()
{

}

void fitFunctionClass::newSplineFunction()
{



 if (parameters[lanSig]>100*parameters[gausSig])
 {
   CreateSplineLandau();
 }
 else if (parameters[lanSig] * 100 < parameters[gausSig]){
   CreateSplineGaus();
 }
 else{
   CreateSplineLandauGauss();
 }






 // axis cumsumLandauGauss_y(convLandau_Gaus_y.size());

}

void fitFunctionClass::newSpline(axis x_axis, axis y_axis)
{
  auto cums = cumsum(y_axis, 1, -1);
  //printVector(cums);
  delete spl;
  spl = new TSpline3("Ilandau", &x_axis[0], &cums[0], cums.size(), "b2e2", 0, 0);
  min_x = x_axis.front();
  max_x = x_axis.back();
}

Double_t fitFunctionClass::operator()(Double_t x)
{

  if (x<min_x)
  {
    return parameters[Ampl];
  }
  else if (x>max_x)
  {
    return 0;
  }
  return spl->Eval(x)*parameters[Ampl];
}

void fitFunctionClass::setNewParameter(Double_t* newParameter)
{
  bool isEqual = true;
  for (size_t i = 0; i < parSize;++i)
  {
    if (parameters[i]!=newParameter[i])
    {
      isEqual = false;
      break;
    }
  }

  
  if (!isEqual)
  {
 //   cout << "new Parameter" << endl;
    for (size_t i = 0; i < parSize; ++i)
    {
      parameters[i] = newParameter[i];

    }
    newSplineFunction();
  }

}

void fitFunctionClass::CreateSplineGaus()
{
  value_t delta_x = parameters[gausSig] / 10;
  auto gaus_x = makeLine(-5 * parameters[gausSig], delta_x, 5 * parameters[gausSig]);
  auto gauss_y = callFunction(&TMath::Gaus, gaus_x, 0, parameters[gausSig]);

  newSpline(gaus_x, gauss_y);
}



void fitFunctionClass::CreateSplineLandau()
{
  value_t delta_x =parameters[lanSig] / 5;

  auto landau_x = makeLine(parameters[lanMp] - 5 * parameters[lanSig], delta_x, parameters[lanMp] + 100 * parameters[lanSig]);
  auto landau_y = callFunction(&TMath::Landau, landau_x, parameters[lanMp], parameters[lanSig]);
  newSpline(landau_x, landau_y);
}

void fitFunctionClass::CreateSplineLandauGauss()
{
  value_t delta_x = __min(parameters[gausSig] / 10, parameters[lanSig] / 5);

  auto gaus_x = makeLine(-5 * parameters[gausSig], delta_x, 5 * parameters[gausSig]);
  auto landau_x = makeLine(parameters[lanMp] - 5 * parameters[lanSig], delta_x, parameters[lanMp] + 100 * parameters[lanSig]);


  auto gauss_y = callFunction(&TMath::Gaus, gaus_x, 0, parameters[gausSig]);
  auto landau_y = callFunction(&TMath::Landau, landau_x, parameters[lanMp], parameters[lanSig]);

  auto convLandau_Gaus_y = convolut(landau_y, gauss_y);
  value_t sum = 0;
  for (auto& e : convLandau_Gaus_y)
  {
    e *= (delta_x*delta_x);
    sum += e;
  }

  //cout << "sum: " << sum << endl;
  assert(abs(sum - 1) < 0.05);
  // printVector(convLandau_Gaus_y);
  // printVector(landau_y);
  // printVector(gauss_y);

  auto convLandau_Gaus_x = makeLine(landau_x.front() + gaus_x.front() - delta_x, delta_x, landau_x.back() + gaus_x.back() + delta_x);
  assert(convLandau_Gaus_x.size() == convLandau_Gaus_y.size());

  newSpline(convLandau_Gaus_x, convLandau_Gaus_y);

}
