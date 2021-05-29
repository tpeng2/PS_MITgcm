% load grid
dyg=rdmds([gcmpath,'DYG']); drF=rdmds([gcmpath,'DRF']); dxg=rdmds([gcmpath,'DXG']);
xg=rdmds([gcmpath,'XG']); yg=rdmds([gcmpath,'YG']); rF=rdmds([gcmpath,'RF']); 
xc=rdmds([gcmpath,'XC']); yc=rdmds([gcmpath,'YC']); rC=rdmds([gcmpath,'RC']); 

hFacC=rdmds([gcmpath,'hFacC']); 
hFacCnan=hFacC;
hFacCnan(hFacC==0)=NaN;
hFacCmsk=hFacCnan;
hFacCmsk(hFacC>0)=1;
Nx=size(dxg,1); Ny=size(dyg,2); Nz=size(drF,3);
dx=Lx/Nx; dy=Ly/Ny;
x=dx/2:dx:dx*Nx-dx/2; y=dy/2:dy:dy*Ny-dy/2; 
%2D grids
[X Y]=meshgrid(x,y); X=X';Y=Y';
xvort=dx:dx:dx*Nx-1; yvort=dy:dy:dy*Ny-1;
% 3D grid spacing
drF_3d = repmat(drF,[Nx,Ny]); 
