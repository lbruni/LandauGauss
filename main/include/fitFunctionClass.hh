#ifndef fitFunctionClass_h__
#define fitFunctionClass_h__

#include "Rtypes.h"
#include <vector>
#define lanMp   0
#define lanSig  1
#define Ampl    2
#define gausSig 3
#define parSize 4
#ifdef WIN32

#define DllExport   __declspec( dllexport )
#else 
#define DllExport   
#endif // WIN32

class TSpline3;
typedef Double_t value_t;
typedef std::vector < value_t > axis;

class DllExport fitFunctionClass{
public:
  fitFunctionClass(){
  
    parameters[gausSig] = 1;
    parameters[lanMp] = 0;
    parameters[lanSig] = 1;
    parameters[Ampl] = 1;
    newSplineFunction();
  }
  void newFunction();
  void newSplineFunction();
  void setNewParameter(Double_t* newParameter);
  Double_t operator()(Double_t x);
private:
  void CreateSplineGaus();
  void CreateSplineLandau();
  void CreateSplineLandauGauss();

  void newSpline(axis x_axis, axis y_axis);
  
  Double_t parameters[parSize];

  TSpline3* spl=nullptr;
  Double_t min_x,max_x;
  
};

Double_t NewLandauGausInt(Double_t *x, Double_t *par);

void SetNewLandauGauss_Setparation(Double_t Separation);



#endif // fitFunctionClass_h__
