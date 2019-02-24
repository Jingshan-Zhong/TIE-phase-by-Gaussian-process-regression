% POISSONFFT  - solves poisson eqn by FFT method 
% based on Frankot and Chellappa algorithm
% modified by Laura Waller, 2009, MIT lwaller@alum.mit.edu
%
% function z = poissonFFT(dzdxy,zpad,regparam)
% inputs: dzdxy = Laplacian of field to be recovered
%         zpad = size of matrix to use in FFT, for zero padding
%         regparam = regularization parameter 
% output: z = solution to poisson equation - assumes periodic boundary
%         conditions

%last updated by Jingshan, Nov 24, 2013

function z = poissonFFT(dzdxy,regparam)
[nx,ny]=size(dzdxy);
rows=nx+100; cols=ny+100;
%[rows,cols]=size(dzdxy);
[wx, wy] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
    ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));
wx = ifftshift(wx); wy = ifftshift(wy);
DZDXY = fft2(dzdxy,rows,cols);
%divide by Poisson soln in Fourier domain
Z = (DZDXY)./(4*pi*pi*(wx.^2 + wy.^2 + regparam));
%Z = (DZDXY).*(4*pi*pi*(wx.^2 + wy.^2))./( (4*pi*pi*(wx.^2 + wy.^2)).^2+ regparam);
z = real(ifft2(Z));    %solution
z=z(1:nx,1:ny);        %get rid of zero padded area
