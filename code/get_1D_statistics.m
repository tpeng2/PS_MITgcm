        % 1D files
        % KE by layer
        KE_3D=1/2.*(UV(:,:,:,1).*UV(:,:,:,1)+UV(:,:,:,2).*UV(:,:,:,2));
        KE_1d=squeeze(nanmean(KE_3D,[1,2]));
        KE1D_xymean(:,i)=KE_1d;
        
        % mean velocity and temperature
        T_mean=squeeze(nanmean(WT(:,:,:,2),[1,2]));
        T_mean_3D=permute(repmat(T_mean,[1,Nx,Ny]),[2,3,1]);
        W_mean=squeeze(nanmean(WT(:,:,:,1),[1,2]));
        W_mean_3D=permute(repmat(W_mean,[1,Nx,Ny]),[2,3,1]);
        U_mean=squeeze(nanmean(UV(:,:,:,1),[1,2]));
        U_mean_3D=permute(repmat(U_mean,[1,Nx,Ny]),[2,3,1]);
        V_mean=squeeze(nanmean(UV(:,:,:,2),[1,2]));
        V_mean_3D=permute(repmat(V_mean,[1,Nx,Ny]),[2,3,1]);

        % perturbation
        U_turb=squeeze(nanmean(UV(:,:,:,1)-U_mean_3D,[1,2]));
        V_turb=squeeze(nanmean(UV(:,:,:,2)-V_mean_3D,[1,2]));
        W_turb=squeeze(nanmean(WT(:,:,:,1)-W_mean_3D,[1,2]));
        T_turb=squeeze(nanmean(WT(:,:,:,2)-T_mean_3D,[1,2]));
    
        % TKE
        U_var2=squeeze(nanmean((UV(:,:,:,1)-U_mean_3D).^2,[1,2]));
        V_var2=squeeze(nanmean((UV(:,:,:,2)-V_mean_3D).^2,[1,2]));
        W_var2=squeeze(nanmean((WT(:,:,:,1)-W_mean_3D).^2,[1,2]));
        T_var2=squeeze(nanmean((WT(:,:,:,2)-T_mean_3D).^2,[1,2]));
        
        
        mat_1D_save=[U_mean,V_mean,W_mean,T_mean,U_turb,V_turb,W_turb,T_turb,U_var2,V_var2,W_var2,T_var2,KE_1d];
        
        fid_stat1d=fopen([path_save,'Stats_1D.',num2str(fileind,'%010d')],'w');
        fwrite(fid_stat1d,mat_1D_save,'float');
        fclose(fid_stat1d);
        
        
