%Inverse of the coeff

function FrequencyMesh=CalFrequency(Img,lambda, ps,dz)
[nx,ny]=size(Img);

dfx = 1/nx/ps;
dfy = 1/ny/ps;
[Kxdown,Kydown] = ndgrid(double(-nx/2:nx/2-1),double(-ny/2:ny/2-1));
k0 = find(Kxdown==0&Kydown==0);
Kxdown = Kxdown*dfx;
Kydown = Kydown*dfy;

FrequencyMesh=lambda*pi*(abs(Kxdown.^2)+abs(Kydown.^2));
FrequencyMesh=FrequencyMesh*dz/(2*pi); % Normalized for sampling step and GP regression 

%
%% show Frequency
%
% DiagF=diag(FrequencyMesh);
% figure(1);
% subplot(2,2,1);
% imagesc(FrequencyMesh);
% axis image;axis off;colormap gray
% title('FrequencyMesh');colorbar
% 
% subplot(2,2,2)
% plot([1:length(DiagF)],DiagF);
% title('Cross line');



