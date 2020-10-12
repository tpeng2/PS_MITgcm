%% Find DiagUV File Names
imageNames_UV = dir(fullfile(gcmpath,['Diag_snaps_UV*.data']));
imageNames_UV = {imageNames_UV.name}';
total_slices_UV=length(imageNames_UV);
for i=1:1:total_slices_UV
imageid_UV(i)=str2num(imageNames_UV{i}(15:end-5));
end

% trim file indices that are out of range
leftid_UV=find(imageid_UV<=indmax_UV);
rightid_UV=find(imageid_UV>=indmin_UV);
trimfid_UV=intersect(leftid_UV,rightid_UV);
if(max(imageid_UV)<indmin_UV || min(imageid_UV)>indmax_UV)
disp('The file indices are not valid.');
end
% define reading length
total_rdslice_UV=length(trimfid_UV);
fileind_UV=imageid_UV(trimfid_UV);
% start and end indices for saving files
startind_UV=fileind_UV(1); endind_UV=fileind_UV(end);
total_day=total_slices_UV*Diag_UV_freq/86400;
startyear=startind_UV*dt/Diag_UV_freq/86400/360;
total_year_UV=total_day/360;
yearind_UV=linspace(startyear,startyear+total_year_UV,total_rdslice_UV);

% interpolation x-z plane
Xg2d=reshape(xg_3d(:,ny/2,:),[nx,nz]); 
ZC2d=reshape(rC_3d(:,ny/2,:),[nx,nz]); 
[Xg2d,ZC2d]=meshgrid(xg(:,1,1),rC);
[Xgint,ZCint]=meshgrid(xg(:,1,1),rC_unif);
dz_unif=mean(diff(rC_unif));
%% Initialization
% meridional transport
Transp=zeros(total_rdslice_UV,1);
Vort_top=zeros(nx,ny,total_rdslice_UV);
Vort_basin=zeros(nx,ny,total_rdslice_UV);
Vort_bot=zeros(nx,ny,total_rdslice_UV);
Vort_xz=zeros(nx,nz_interp,total_rdslice_UV);
inctime=0;

%%
dind=0;ind_old=0;
mkdir(path_output,'images')
for i=1:1:total_rdslice_UV
    ind=fileind_UV(i);
    disp(['processing UV',num2str(i/total_rdslice_UV*100,'%.4f'),'%, f_ind:',num2str(ind),' / ',num2str(fileind_UV(end))]);
    if(i>1) dind=ind-ind_old;end
    % load temporary U
    UV=rdmds([gcmpath,'Diag_snaps_UV.',num2str(ind,'%010d')]);
    U=UV(:,:,:,1); V=UV(:,:,:,2);
    % calculate zonal transport
    Transp(i)=mean(sum(sum((U.*dyg_3d).*drF_3d)));
    % calculate relative vorticity (k-component)
    Utop=U(:,:,1); Vtop=V(:,:,1);
    Vort_top(:,:,i)=curl(Utop,Vtop)/dx;
    Ubasin=U(:,:,72); Vbasin=V(:,:,72);
    Vort_basin(:,:,i)=curl(Ubasin,Vbasin)/dx;
    Ubot=U(:,:,end-10); Vbot=V(:,:,end-10);
    Vort_bot(:,:,i)=curl(Ubot,Vbot)/dx;
    % calculate j-componenet vorticity (-dw/dx+du/dz)
    % at (:,ny/2,:) plane
     Uctr_xz=reshape(U(:,ny/2,:).*Zwetdry(:,ny/2,:),[nx,nz]); 
     Wctr_xz=reshape(W(:,ny/2,:).*Zwetdry(:,ny/2,:),[nx,nz]); 
    % interpolate 
     Uctr_xz_unif=interp2(Xg2d,ZC2d,Uctr_xz',Xgint,ZCint);
     Wctr_xz_unif=interp2(Xg2d,ZC2d,Wctr_xz',Xgint,ZCint);
     Vort_xz(:,:,i)=Vorticity_AnyStencil(dx,dz_unif,Uctr_xz_unif,Wctr_xz_unif,3)';
    if(plot_vort==1) 
     %plotting slides
        f1a=figure('visible','off');
        set(gcf,'Units','pixels','Position',[50 50 960 640]);
        subplot(1,2,1)
        h=pcolor(X./1e3,Y./1e3,Vort_top(:,:,i)./f_2d);colorbar,caxis([-4e-1,4e-1]);
        set(h, 'EdgeColor', 'none');
        xlabel('x [km]');xlabel('y [km]'); title('\zeta/f at the surface');
        %annotation
        annotation('textbox',[0.02,0.05,0.10,0.06],...
            'string',{['t_0+',num2str(inctime,'%.3f'),'day']}); axis off
        subplot(1,2,2)
        h=pcolor(X./1e3,Y./1e3,Vort_bot(:,:,i)./f_2d);colorbar,caxis([-1e-1,1e-1]);
        set(h, 'EdgeColor', 'none');
        xlabel('x [km]');xlabel('y [km]'); title('\zeta/f near topography')
        saveas(f1a,[path_output,'./images/','vort_DiagUV_',num2str(ind),'.jpg'])
        close(f1a)
    end %if plot_vort==1
    %time interval
    inctime=dind*dt/86400*i;
    ind_old=ind;
end
%%
save([path_output,'VORT_',num2str(startind_UV,'%010d'),'-',num2str(Diag_UV_freq/86400,'%010d'),'d-',num2str(endind_WT,'%010d'),'.mat'],...
    'Vort_top','Vort_bot','Vort_xz','Transp','Xgint','ZCint','X','Y');

%%
close all
f1=figure(3);
plot(yearind_UV,Transp/1e6,'-.'),hold on
plot(yearind_UV,ones(length(Transp),1)*mean(Transp'/1e6),'--')
legend('Transport',['mean=',num2str(mean(Transp/1e6)),'Sv'])
title('Zonal transport')
ylabel('Volume flux [Sv]')
xlabel('Years')
saveas(f1,[path_output,'Zonal_transp_',num2str(total_year_UV,'%.3f'),'.png'])
