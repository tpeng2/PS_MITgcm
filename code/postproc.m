clear,clc
set(0,'defaultaxesfontsize',14)
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultlinelinewidth',1.2) 
%%
plot_vort=0;
%% Time stepping
dt=120;% per step
indmin_UV=0; indmax_UV=02700000;
indmin_WT=0; indmax_WT=02700000;
Diag_WT_freq=5184000;
Diag_UV_freq=5184000;

%% model parameters
f0=-1e-4; beta=1.5E-11; 
Lx=480e3; Ly=960e3; Lz=2e3;
%% load z and calculate topography
addpath('../../code/')
groupname='Top_CH-ACC';
casename='woce2km_288c';
case_path='/scratch/tpeng/MITgcm_cases/';
post_result_path='/scratch/tpeng/MITgcm_post/results/';
post_code_path='/scratch/tpeng/MITgcm_post/code/';
gcmpath=[case_path,casename,'/'] 
path_output=[post_result_path,'/',groupname,'/',casename,'/'] 
mkdir(path_output);
syscommand=['sh ',post_code_path,'/sh_scrpt/copy_grids.sh ',groupname,' ',casename,' ',case_path,' ',post_result_path];
system(syscommand)
%%
% construct grids
dyg=rdmds([gcmpath,'DYG']); drF=rdmds([gcmpath,'DRF']); dxg=rdmds([gcmpath,'DXG']);
xg=rdmds([gcmpath,'XG']); yg=rdmds([gcmpath,'YG']); rF=rdmds([gcmpath,'RF']); 
xc=rdmds([gcmpath,'XC']); yc=rdmds([gcmpath,'YC']); rC=rdmds([gcmpath,'RC']); 

hFacC=rdmds([gcmpath,'hFacC']); 
hFacCnan=hFacC;
hFacCnan(hFacC==0)=NaN;
hFacCmsk=hFacCnan;
hFacCmsk(hFacC>0)=1;
nx=size(dxg,1); ny=size(dyg,2); nz=size(drF,3);
dx=Lx/nx; dy=Ly/ny;
x=dx/2:dx:dx*nx-dx/2; y=dy/2:dy:dy*ny-dy/2; 
%2D grids
[X Y]=meshgrid(x,y); X=X';Y=Y';
xvort=dx:dx:dx*nx-1; yvort=dy:dy:dy*ny-1;
% 3D grid spacing
dxg_3d = repmat(dxg,[1,1,nz]);xg_3d = repmat(xg,[1,1,nz]);
dyg_3d = repmat(dyg,[1,1,nz]);yg_3d = repmat(yg,[1,1,nz]);
drF_3d = repmat(drF,[nx,ny]); rF_3d = repmat(rF,[nx,ny]);
rC_3d = repmat(rC,[nx,ny]);
drF_3dFac=drF_3d.*hFacC; 
vol_grid=(dxg_3d.*dyg_3d).*drF_3d;
vol_gridFac=(dxg_3d.*dyg_3d).*drF_3d.*hFacCnan;
% Uniform z-grid for interpolation
nz_interp=2*nz;
rC_unif=linspace(0,-4000,nz_interp); 
rC_unif=reshape(rC_unif,[1,1,nz_interp]);
rC_3d_unif=repmat(rC_unif,[nx,ny]);

% Total volume and horizontal area
Vol=(nx*dx*ny*dy*Lz);horArea=nx*dx*ny*dy;
% vertical grids
z_cntr=1/2*drF(1)+cumsum(-drF); 
z_cntr=squeeze(z_cntr);% horizontal domain
z_upper=zeros(size(drF,3),1);
z_upper(2:end)=cumsum(-drF(1:end-1)); z_upper=squeeze(z_upper);% horizontal domain
z_lower=zeros(size(drF,3),1); z_lower(end)=-4000;z_lower(1:end-1)=cumsum(-drF(1:end-1));
z_lower=squeeze(z_lower);% horizontal domain
zFac=cumsum(drF_3dFac,3); %effective z
% coriolis map
f=f0-y.*beta;f_2d=repmat(f,[nx,1]);
%% load topography
ieee='b';
accuracy='real*8';
fid_top=fopen([gcmpath,'topog_',num2str(nx),'_',num2str(ny),'_sigma_300.box'],'r',ieee) 
Zb=fread(fid_top,nx*ny,'real*8');
Zb=reshape(Zb,[nx,ny]);
Zwetdry=zeros(nx,ny,nz);Zwetdry(Zwetdry==0)=NaN;
%% find last full wet grid index and the wet fraction below it
iZ2d_up=zeros(size(Zb)); %1st grid higher than the boudnary
Zfrac2d_wet=zeros(size(Zb)); %wet volume
% find bottom 
for i=1:nx;
    for j=1:ny;
        izbot= min(find(Zb(i,j)>=z_lower));
        iZ2d_up(i,j)=izbot-1;
        dz_upper=z_upper(izbot)-Zb(i,j);
        fracz_wet=dz_upper/drF(izbot);
        if fracz_wet<0.2
            fracz_wet=0;    
        end
        Zfrac2d_wet(i,j)=fracz_wet;
        % weighting grid for full wet cells
        Zwetdry(i,j,1:iZ2d_up(i,j))=1;
        % weighting grid for fractional wet cell
        Zwetdry(i,j,izbot)=fracz_wet;
    end
end
iZ2d_frac=iZ2d_up-1;
ind_dry=find(Zfrac2d_wet==0);
iZ2d_frac(ind_dry)=iZ2d_up(ind_dry);


%% 
read_DiagWT
%%
read_DiagUV
%%
if plot_vort==1
    vorticity_movie
end
