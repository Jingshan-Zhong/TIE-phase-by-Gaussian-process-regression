****************************** INTRODUCTION**************************

Note any commercial use or modification of the code is not allowed without the permission of the authors.
In the GP TIE algorithm, we perform Gaussian process regression over the defocused intensity images [in Frequency domain] to estimate the intensity axial derivative, which is used to recover phase by the transport intensity of equation (TIE) method. GP TIE alleviates the nonlinearity error in the derivative estimation by using the prior knowledge of how intensity varies with defocus propagation in the spatial frequency domain. It doesn’t require the intensity images to be equally spaced, so the input intensity stack can be exponentially spaced, which is shown to be an efficient scheme to transfer the phase information into the measured intensity. For more details, please see the reference paper.

************************HOW TO USE THE CODE**************************
How to run on the example data set:
1)	Open Main_GPTIE.m and run in Matlab. The example data set will be automatically loaded.
How to run on your own data set:
1)	Prepare your data set followed the format of example data set ‘SampleData2.mat’. Make sure that the variable names are same as the example data and the unit of the variable is meter.
2)	Open Main_GPTIE.m.
3)	Load in your data and run.
4)	Tune Poisson solver regularization parameter (regparam) if necessary.

The Input parameters:
Inten: defocused intensity stack measured along the propagation axis
zvec: the positions of the measured intensity images [m]
ps: pixel size [m]
lambda: wavelength [m]
 
Output:
RePhase1: Recovered phase [radian]
 
Default parameters:
zfocus1: the position of the focal plane, which is defaulted as zfocus1=0 [m].
regparam=5*10^-6; Poisson solver regularization for GP TIE

Note: 
1) The variable zvec needs to be a COLUMN vector.
2) The current version only works for pure phase recovery, which means the intensity
at focus is assumed to be constant

***************************REFERENCE PAPER ***************************
Zhong Jingshan, Rene A. Claus, Justin Dauwels, Lei Tian, and Laura Waller, "Transport of Intensity phase imaging by intensity spectrum fitting of exponentially spaced defocus planes," Opt. Express 22, 10661-10674 (2014)  



