%% Find DiagWT File Names
imageNames_WT = dir(fullfile(gcmpath,['Diag_snaps_WT*.data']));
imageNames_WT = {imageNames_WT.name}';
total_slices_WT=length(imageNames_WT);
for i=1:1:total_slices_WT
imageid_WT(i)=str2num(imageNames_WT{i}(15:end-5));
end

% trim file indices that are out of range
leftid_WT=find(imageid_WT<=indmax_WT);
rightid_WT=find(imageid_WT>=indmin_WT);
trimfid_WT=intersect(leftid_WT,rightid_WT);
if(max(imageid_WT)<indmin_WT || min(imageid_WT)>indmax_WT)
disp('The file indices are not valid.');
end
% define reading length
total_rdslice_WT=length(trimfid_WT);
fileind_WT=imageid_WT(trimfid_WT);
% start and end indices for saving files
startind_WT=fileind_WT(1); endind_WT=fileind_WT(end);
total_day=total_slices_WT*Diag_WT_freq/86400;

%% Read T from Diagnostics
Tsouth_avg=zeros(nz,total_rdslice_WT);
Tcentr_avg=zeros(nz,total_rdslice_WT);
Tnorth_avg=zeros(nz,total_rdslice_WT);
Tz1d_avg=zeros(nz,total_rdslice_WT);
Tvolmean=zeros(total_rdslice_WT,1);
Ld_mean=zeros(total_rdslice_WT,1);
dTdzm_south=zeros(nz-1,total_rdslice_WT);
dTdzm_north=zeros(nz-1,total_rdslice_WT);
dTdzm_centr=zeros(nz-1,total_rdslice_WT);
dTdzm_vol=zeros(nz-1,total_rdslice_WT);
f0=-7.2722e-5; beta=1.5E-11; f=f0-y.*beta;
f_2d=repmat(f,[nx,1]);
%%
%startind_WT=min(find(ind>=indmin_WT));
disp(['Total slides to read: ', num2str(total_rdslice_WT,'%06d')])
for i=1:1:total_rdslice_WT
        ind=fileind_WT(i);
        disp(['Processing WT: ',num2str(i/total_rdslice_WT*100,'%.3f'),'%, f_ind:',num2str(ind),' / ',num2str(fileind_WT(end))]);
        % for file name identifier
        if(i==1), startind_WT=ind;elseif (i==total_slices_WT), endind_WT=ind;end
%         disp(['No. ',num2str(i,'%06d'),', file_ind: ',num2str(ind,'%010d')])
        % load temporary U
        WT=rdmds([gcmpath,'Diag_snaps_WT.',num2str(ind,'%010d')]);
        W=WT(:,:,:,1); T=WT(:,:,:,2);
        T(T==0)=NaN;
        % volume averaged temperature
        Tvol=T.*vol_gridFac;
        Tzmean=nanmean(nanmean(T.*hFacCmsk,1));Tzmean=squeeze(Tzmean);
        Tvolmean(i)=nansum(nansum(nansum(Tvol)))/Vol;
        Ld_mean(i)=get_Lrho_T(Tzmean,squeeze(drF),f0);
        % calculate vorticity (surface layer)
        wallind_S=3; wallind_N=ny-2; wallind_C=ny/2;
        Tsouth=T(:,wallind_S,:); Tnorth=T(:,wallind_N,:); Tcentr=T(:,wallind_C,:);
        Tsouth=squeeze(Tsouth);
        Tcentr=squeeze(Tcentr);
        Tnorth=squeeze(Tnorth);    
        Tsouth_avg(:,i)=nanmean(Tsouth);
        Tcentr_avg(:,i)=nanmean(Tcentr);
        Tz1d_avg(:,i)=nanmean(nanmean(T));
        dzF=diff(zFac,[],3);dzF(dzF==0)=NaN;
        Tnorth_avg(:,i)=nanmean(Tnorth);
        zF_S=squeeze(zFac(:,wallind_S,:));   dzF_south=squeeze(dzF(:,wallind_S,:));
        zF_C=squeeze(zFac(:,wallind_C,:));   dzF_centr=squeeze(dzF(:,wallind_S,:));
        zF_N=squeeze(zFac(:,wallind_N,:));   dzF_north=squeeze(dzF(:,wallind_S,:));
        dTdzm_south(:,i)= nanmean(diff(Tsouth,[],2)./dzF_south);
        dTdzm_centr(:,i)= nanmean(diff(Tcentr,[],2)./dzF_centr);
        dTdzm_north(:,i)= nanmean(diff(Tnorth,[],2)./dzF_north);
        difTz=diff(T,[],3)./dzF;
        dTdzm_vol(:,i)=nanmean(nanmean(difTz,1),2);
end

%%
save([path_output,'T_NS_',num2str(startind_WT,'%010d'),'-',num2str(Diag_WT_freq/86400,'%010d'),'d-',num2str(endind_WT,'%010d'),'.mat'],'Tz1d_avg','Tsouth_avg','Tcentr_avg','Tnorth_avg','Tvolmean','startind_WT','Diag_WT_freq','endind_WT','dt','drF','dxg','dyg','dTdzm_south','dTdzm_centr','dTdzm_north','dTdzm_vol','Ld_mean','fileind_WT')
