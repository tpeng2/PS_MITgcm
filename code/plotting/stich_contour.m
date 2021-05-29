matpath='/home/tpeng/MITgcm/postproc/results/Top_CH-ACC/run4km_160c';
%%
nz=100;
len1=395;
len2=0;
len3=0;
len=len1+len2+len3;
%create a long array for time
FID=zeros(1,len);
TCNTR=zeros(nz,len);
LDMEAN=zeros(1,len);
TVOLM=zeros(1,len);
DTDZMVOL=zeros(nz-1,len);
YEAR=zeros(1,len);
%% fill ind
fillind=1:len1;
% fillind=len1+1:len1+len2;
% fillind=len1+len2+1:len1+len2+len3;
indinc_WT=5040; Diag_WT=89750;
% yearind_WT=fileind_WT*dt/86400/360;
YEAR=linspace(2,len*dt*indinc_WT/86400/360,len);
%%
FID(fillind)=fileind_WT;
TCNTR(:,fillind)=Tcentr_avg;
LDMEAN(fillind)=Ld_mean;
TVOLM(fillind)=Tvolmean;
DTDZMVOL(:,fillind)=dTdzm_vol;
%%
[YR2d,Z2d]=meshgrid(YEAR,z_cntr);
[dzYR2d Zdz2d]=meshgrid(YEAR,z_dz);
%% 1-D plot
dTdzm1d=nanmean(dTdzm_vol,2);
Tz1dm=nanmean(Tz1d_avg,2);
line(Tz1dm,z_cntr,'Color','r')
ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(dTdzm1d,z_dz,'Parent',ax2,'Color','k')

%%
figure(1)
set(0,'defaultfigurecolor','w')
nlevel=25;
lv_min=1e-4; lv_max=3e-2;
neg_levels=-logspace(log10(lv_min),log10(lv_max),nlevel);
[C,h]=contourf(dzYR2d,Zdz2d,dTdzm_vol,neg_levels,'ShowText','off')
cb=colorbar
set(cb,'Ticks')
xlabel('Year')
ylabel('z[m]')
title('Horizontally-temporarily averaged <dT/dz> [\circ C m^{-1}]')
%%
figure(2)
set(0,'defaultfigurecolor','w')
nlevel=30;
Tmin=2e-2;
Tmax=max(max(Tcentr_avg))+0.1;
levels=logspace(log10(Tmin),log10(Tmax),nlevel);

[C,h]=contourf(YR2d,Z2d,Tz1d_avg,levels,'ShowText','off')
cb=colorbar
% set(cb,'Ticks')
xlabel('Year')
ylabel('z[m]')
title('Horizontally-temporarily averaged <T>(z) [\circ C]')
%%
figure(2)
set(0,'defaultfigurecolor','w')
nlevel=30;
Tmin=1e-4;
Tmax=max(max(Tcentr_avg))+0.1;
levels=logspace(log10(Tmin),log10(Tmax),nlevel);

[C,h]=contourf(YR2d,Z2d,Tz1d_avg,levels,'ShowText','off')
cb=colorbar
% set(cb,'Ticks')
xlabel('Year')
ylabel('z[m]')
title('<T>(z) @center x-z plane [\circ C]')
