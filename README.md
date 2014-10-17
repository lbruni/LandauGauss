LandauGauss
===========

Script that fits a given histogram with a Landau Gauss Function


# Installation
To install the program you need a c++ 11 compiler, Cmake and CERN ROOT.

After downloading create with cmake the make/Project files for your build system. 
> git clone <url>  <your_path>\landau_gauss
> cd <your_path>\landau_gauss\build
> cmake ..

On Linux:
> make  
> make install
On Windows:
> msbuild INSTALL.vcxproj /P:configuration=release

# Known Problems 
On Linux cmake has some problems finding ROOT. I don’t know  yet how to solve it. 

# Using

Example code: 

>  landgausFit fit;
  
>  fit.setStartAmplitude(100);
>  fit.setLimits_Amplitude(95, 100);
>  fit.setLimits_GaussSigma(0, 200);
>  fit.setLimits_LandauMP(0, 500);
>  fit.setLimits_LandauSigma(0, 200);
>  fit.setStartGaussSigma(27);
>  fit.setStartLandauMean(100);
>  fit.setStartLandauSigma(7);
>  fit.setFitRange(0,250);
>  TGraph * g; // create a data set
>  fit(g);  // fit the data set

>  fit.printResults(cout);  // you can also print the result to file with this function 
>  fit.DrawfitFunction(); // draws the dataset with the fit function 
>  fit.saveFitToFile(“test.txt”) // saves details about the fit

