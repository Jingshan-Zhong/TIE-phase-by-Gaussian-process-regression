
function [dIdz Coeff Coeff2]=GPRegression(Ividmeas, zfocus,z,Sigmaf,Sigmal,Sigman)

%Gaussion processs regression

[Nx,Ny,Nz]=size(Ividmeas);
%z=([1:Nz]');
VectorOne=ones(Nz,1);
KZ=VectorOne*z'-z*VectorOne';

K=Sigmaf*(exp(-1/2/Sigmal*KZ.^2));
L=chol(K+Sigman*eye(Nz));

z2=zfocus;

Nz2=length(z2);
VectorOne2=ones(Nz2,1);
KZ2=VectorOne*z2'-z*VectorOne2';

% for first derivative
D=Sigmaf*(exp(-1/2/Sigmal*(KZ2).^2))/(-Sigmal).*(KZ2);
%Coeff=D'/(K+Sigman*eye(Nz));
Coeff=D'/L/L';% use /L/L' to be more stable to the matrix inversion

% for regression
D2=Sigmaf*(exp(-1/2/Sigmal*(KZ2).^2));
%Coeff2=D2'/(K+Sigman*eye(Nz));
Coeff2=D2'/L/L';

%Coeff=Coeff/sum(Coeff2);
%Coeff2=Coeff2/sum(Coeff2); %Normalized,  Necessary More consideration later.

dIdz=zeros(Nx,Ny);

for k=1:Nz
    dIdz=dIdz+Ividmeas(:,:,k)*Coeff(k);
end


