% search file name
fname_UV_all = dir(fullfile(path_diag,['Diag_snaps_UV*.data']));
fname_UV_all = {fname_UV_all.name};
fname_WT_all = dir(fullfile(path_diag,['Diag_snaps_WT*.data']));
fname_WT_all = {fname_WT_all.name};
fname_PH_all = dir(fullfile(path_diag,['PHIHYD*.data']));
fname_PH_all = {fname_WT_all.name};
% %   make it cell
% fname_UV_all={fname_UV_all}';
% fname_WT_all={fname_WT_all}';
%   total slide UV
total_slices_UV=length(fname_UV_all);
total_slices_WT=length(fname_WT_all);
total_slices_PH=length(fname_PH_all);
imageid_UV=zeros(size(fname_UV_all));
imageid_WT=zeros(size(fname_WT_all));
imageid_PH=zeros(size(fname_PH_all));
%   get file indices
for i=1:1:total_slices_UV
        disp(fname_UV_all{i}(15:end-5))
        imageid_UV(i)=str2num(fname_UV_all{i}(15:end-5));
end
for i=1:1:total_slices_WT
        imageid_WT(i)=str2num(fname_WT_all{i}(15:end-5));
        imageid_PH(i)=str2num(fname_PH_all{i}(8:end-5));
end

if ~exist(path_save, 'dir')
           mkdir(path_save);
end

%   trim down file indices
%   UV
leftidind_UV=find(imageid_UV>=startind);
rightidind_UV=find(imageid_UV<=endind);
trimfid_UV=intersect(imageid_UV(leftidind_UV),imageid_UV(rightidind_UV));
if(length(trimfid_UV)<1)
        disp('The file indices are not valid. UV');
end
%   WT
leftidind_WT=find(imageid_WT>=startind);
rightidind_WT=find(imageid_WT<=endind);
trimfid_WT=intersect(imageid_WT(leftidind_WT),imageid_WT(rightidind_WT));
if(length(trimfid_WT)<1)
        disp('The file indices are not valid. WT');
end

%%
slides_ind=trimfid_UV;
Nfiles=length(trimfid_UV);

KE1D_xymean=zeros(Nz,Nfiles);
Transp=zeros(1,Nfiles);

for i=1:Nfiles
        fileind=slides_ind(i);
        tic
        disp(['Reading #',num2str(fileind,'%010d')])
        UV=rdmds([path_diag,'Diag_snaps_UV.',num2str(fileind,'%010d')]);
        WT=rdmds([path_diag,'Diag_snaps_WT.',num2str(fileind,'%010d')]);
        PH=rdmds([path_diag,'PHIHYD.',num2str(fileind,'%010d')]);

        x_vec=1:Nx;
        y_vec=1:Ny;
        
        disp('Read 2D slides')
        % surface slides
        Usurf_mat=squeeze(UV(:,:,1,1));
        read_save_binary_2D(Usurf_mat,'Usurf',path_save,fileind);

        Vsurf_mat=squeeze(UV(:,:,1,2));
        read_save_binary_2D(Vsurf_mat,'Vsurf',path_save,fileind);
        
        Tsurf_mat=squeeze(WT(:,:,1,2));
        read_save_binary_2D(Tsurf_mat,'Tsurf',path_save,fileind);

        Psurf_mat=squeeze(PH(:,:,1));
        read_save_binary_2D(Psurf_mat,'Psurf',path_save,fileind);

        % sub-ekman slides
        Usekm_mat=squeeze(UV(:,:,5,1));
        read_save_binary_2D(Usekm_mat,'Usekm',path_save,fileind);

        Vsekm_mat=squeeze(UV(:,:,5,2));
        read_save_binary_2D(Vsekm_mat,'Vsekm',path_save,fileind);
        
        Tsekm_mat=squeeze(WT(:,:,5,2));
        read_save_binary_2D(Tsekm_mat,'Tsekm',path_save,fileind);

        Psekm_mat=squeeze(PH(:,:,5));
        read_save_binary_2D(Psekm_mat,'Psekm',path_save,fileind);

        % sub-surface w-velocity
	    Wssurf_mat=squeeze(WT(:,:,2,1));
        read_save_binary_2D(Wssurf_mat,'Wssurf',path_save,fileind);
        
        % central depth
        Ucnt_mat=squeeze(UV(:,:,60,1));
        read_save_binary_2D(Ucnt_mat,'Ucnt',path_save,fileind);
        
        Vcnt_mat=squeeze(UV(:,:,60,2));
        read_save_binary_2D(Vcnt_mat,'Vcnt',path_save,fileind);

        Tcnt_mat=squeeze(WT(:,:,60,2));
        read_save_binary_2D(Tcnt_mat,'Tcnt',path_save,fileind);

        Pcnt_mat=squeeze(PH(:,:,60));
        read_save_binary_2D(Pcnt_mat,'Pcnt',path_save,fileind);

        Wcnt_mat=squeeze(WT(:,:,60,1));
        read_save_binary_2D(Wcnt_mat,'Wcnt',path_save,fileind);
        
        % central-channel (y) x-z field
        Uxzcnt_mat=squeeze(UV(:,floor(Ny/2),:,1));
        read_save_binary_2D(Uxzcnt_mat,'Uxzcnt',path_save,fileind);

        Vxzcnt_mat=squeeze(UV(:,floor(Ny/2),:,2));
        read_save_binary_2D(Uxzcnt_mat,'Vxzcnt',path_save,fileind);

        Wxzcnt_mat=squeeze(WT(:,floor(Ny/2),:,1));
        read_save_binary_2D(Wxzcnt_mat,'Wxzcnt',path_save,fileind);

        Txzcnt_mat=squeeze(WT(:,floor(Ny/2),:,2));
        read_save_binary_2D(Txzcnt_mat,'Txzcnt',path_save,fileind);

        toc
        tic
        % 1D statistics
        disp('Get 1D statistics')
        get_1D_statistics

        % get Transport
%        Udydz=bsxfun(@times,UV(:,:,:,1).*hFacC,dyg).*drF_3d;
%        Transp(i)=nanmean(nansum(nansum(Udydz,2),3));
        disp(['Transport: ',num2str(i),',',num2str(Transp(i))])    
        disp(['U_surf, V_surf, and T_surf are saved for #',num2str(fileind,'%010d')])
        toc
        clear UV, clear WT,
end
tic 
save([path_save,'KE1D_xymean_',num2str(slides_ind(1)),'_',num2str(slides_ind(end)),'.mat'],'KE1D_xymean');
toc

