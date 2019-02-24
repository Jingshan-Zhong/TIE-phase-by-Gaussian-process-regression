function dIdzC=CombinePhase(dIdzStack, Frq_cutoff,FrequencyMesh,CoeffStack,Coeff2Stack)
%Cutoff is the cutoff area in CosH
F = @(x) ifftshift(fft2(fftshift(x)));
Ft = @(x) ifftshift(ifft2(fftshift(x)));

[Nx,Ny,Nsl]=size(dIdzStack);

dIdzC_fft=zeros(Nx,Ny);
Maskf=zeros(Nx,Ny);


%figure;
f0=0;
f1=0;
for k=1:Nsl
    dIdz=dIdzStack(:,:,k);
    dIdz_fft=F(dIdz);
    
    f1=Frq_cutoff(k);
    Maskf=zeros(Nx,Ny);
    Maskf(find(FrequencyMesh<=f1&FrequencyMesh>f0))=1;
    f0=f1;
    dIdzC_fft=dIdzC_fft+(dIdz_fft.*Maskf); %Only update in the update area
     
     %{
%     Coeff=CoeffStack(:,k);
%     Coeff_f=fftshift(fft(Coeff));
%   
%   
%     Coeff2=Coeff2Stack(:,k);
%     Coeff2_f=fftshift(fft(ifftshift(Coeff2)));
    
%     Nz=length(Coeff);
%     subplot(3,2,1);
%     imagesc(Maskf);
%     axis image;axis off;colormap gray
%     title(sprintf('Mask %d',k));colorbar
%     
%     
%     subplot(3,2,3);
%     plot([1:Nz],Coeff);
%     title(sprintf('SumCoeff=%f,SumCoeff^2=%f',sum(Coeff),sum(Coeff.^2)));
%     
%     subplot(3,2,4);
%     plot([1:Nz],abs(Coeff_f));
%     title('Coeff of derivative in Fourier domain');
%     
%     subplot(3,2,5);
%     plot([1:Nz],Coeff2);
%     title(sprintf('SumCoeff=%f,SumCoeff^2=%f',sum(Coeff2),sum(Coeff2.^2)));
%     
%     subplot(3,2,6);
%     plot([1:Nz],abs(Coeff2_f));
%     title('Coeff in Fourier domain');
%     
%      
%     pause(0.05);
      
    %}
end

dIdzC=real(Ft(dIdzC_fft));