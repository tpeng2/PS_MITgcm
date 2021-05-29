% input 3D field
function read_save_binary_2D(M_mat,M_fname,path_save,fileind)
fid_M=fopen([path_save,M_fname,'.',num2str(fileind,'%010d')],'w');
fwrite(fid_M,M_mat,'float');
fclose(fid_M);
return
