%GP TIE code (By Jingshan Zhong)
%note no commercial use or modification of this code is allowed without
%permission of the authors.
%reference:  Transport of Intensity phase imaging by intensity spectrum
%fitting of exponentially spaced defocus planes (Optics Express 2014 accepted)
%contact: www.dauwels.com  or  www.laurawaller.com

%**************************************************************************
% Input paramters:
% Inten: defocused intensity stack measured along the propgagation axis
% zvec: the positions of the images been measured. [meters]
% ps: pixel size [meters]
% lambda: wavelength [meters]

%Output:
%RePhase1: Recovered phase

%Default parameters:
% zfocus1: the position of focal plane is defaulted as zfocus1=0
%regparam=5*10^-6; %Poisson solver regularization for GP TIE
%**************************************************************************

%Important: 
%1) The variable z need to be a COLUMN vector.
%2) This code only works for pure phase recovery, which means the intensity
%3£©The intensity at focus is assumed to be constant

%% loading experimental data
close all;clear all;
load SampleData2;
Ividmeas1=Inten; %Inten stores the intensity images of defocal planes along the direction of light propagation
z1=zvec; %1) z1 is a column vector; 2) positons of the defocal planes in Inten
clear('Inten','zvec');

%default parameters
zfocus1=0; % in the data set, the infocus plane is defaulted at zfocus1=0;
regparam=5*10^-6; %Poisson solver regularization for GP TIE


%% run Gaussian Process
[Nx,Ny,Nz]=size(Ividmeas1);
Nsl=50; %Nsl is defaulted as 50. In order to reduce the computational complexity, we divide the frequency to Nsl bins. 
%For the frequency within one bin, it shares same frequency threshold and same hyper-parameters in GP regression.

RePhase1=RunGaussionProcess(Ividmeas1,zfocus1,z1,lambda,ps,Nsl,regparam);% recovered phase
RePhase1=RePhase1/mean(mean(mean(Ividmeas1))); % normalized phase by mean of intensity.

%% Show result for experimental data
%
figure;
imagesc(RePhase1);
axis image;axis off;colormap gray
title('Recovered Phase by GP TIE');colorbar
%}


%% Show intensity images
%{
figure;
for k=1:length(zvec)
    
    imagesc(Ividmeas1(:,:,k));
    axis image;axis off;colormap gray
    title(sprintf('Noised intensity by this code at z step %d',k));colorbar
   
    pause(0.1)
end
%}
