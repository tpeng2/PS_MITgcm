function [L_d]=get_Lrho_T(T_in,drF,f0)
nz=length(drF);
z_cntr=1/2*drF(1)+cumsum(-drF); z_cntr=reshape(z_cntr,[1,nz]);% horizontal domain
z_lower=cumsum(-drF);
profile_tot=zeros(1,nz);
profile_tot=T_in; 

T_in(T_in==NaN)=0;
% find non-zero profile
profile=T_in(T_in~=0);
nz2=sum(profile_tot~=0);

% extend by one point
profile(nz2+1)=profile(nz2); %nz+1

% Total depth
Htot=z_lower(end);

% Gravity and reduced gravity
g_0=9.81;
g=-diff(profile)*(g_0*2E-4);
h=drF;

% coumpute G
% F=1/(Rossby Radius)^2
F_layer=zeros(nz,nz); % dim: Nz x Nz
% Zeroth layer
F_layer_zero=(f0^(2))/(h(1)*g_0);
% First mode
F_layer(1,1) = (f0^(2))/(h(1)*g(1));
for k=2:nz2    
    %upper diagonal
        F_layer(k,k-1) = (f0^(2))/(h(k)*g(k-1));
    % Diagonal
    F_layer(k,k) = (f0^(2))/(h(k)*g(k));
end
% Last mode forced at 0
F_layer(nz2,nz2) = 0;  

% !F_mode = eigenvalues of "mode" system of equation : q_mode = laplacian(psi_mode)+F*psi_mode
M=zeros(nz2,nz2); 
%First value
M(1,1) = - (F_layer_zero(1)+F_layer(1,1));
M(2,1) =    F_layer(2,1);
for k=2:nz2-1
    M(k-1,k) =    F_layer(k-1,k-1);
    M(k,k)   = - (F_layer(k,k-1) + F_layer(k,k));
    M(k+1,k) =    F_layer(k+1,k);
end 
%Last group
  k = nz2;
M(k-1,k) =   F_layer(k-1,k-1);
M(k,k)   = - F_layer(k,k-1);       

MT=M';

% Replace NaN
MT(isnan(MT))=0;
MT(isinf(MT))=0;
[V,D] = eig(full(MT));
      

D=diag(D);
size_V=size(V);
size_D=size(D);


Rossby=(sqrt(-1./D(:)) / 1000)';


Rossby =sort(Rossby);

L_d=Rossby(nz2-1);
end