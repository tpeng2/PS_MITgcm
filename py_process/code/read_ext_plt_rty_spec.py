
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import matplotlib
matplotlib.use('Agg')
home_dir=os.path.expanduser("~")
package_path=home_dir+'/MITgcm_post/py_process/ppmod/'
os.chdir(package_path)
sys.path.append(package_path)
from PPFCN import *
import PPFCN.ppfcns as ppf
from scipy import interpolate
#%%
import argparse
import shutil


parser = argparse.ArgumentParser(description="Read ExfSuf")
parser.add_argument('groupname', metavar="GROUPNAME", help="Group name (groupname)")
parser.add_argument('casename', metavar="CASENAME", help="Case name (casename)")
parser.add_argument('start_ind', metavar="START_IND", help="Start index (start_ind)")
parser.add_argument('end_ind', metavar="END_IND", help="Start index (end_ind)")
parser.add_argument('loc_str', metavar="LOC_STR", help="Location (e.g.,  'surf', 'sekm', or 'cnt'.")

args = parser.parse_args()

groupname=args.groupname
casename=args.casename
start_ind=int(args.start_ind)
end_ind=int(args.end_ind)
loc_str=args.loc_str
print('Read groupname as '+args.groupname)
print('Read casename as '+args.casename)
print('Read start_ind as '+args.start_ind)
print('Read end_ind as '+args.end_ind)
print('Read loc_str as '+args.loc_str)
# groupname=''
# casename='rel_uvar_9'

# start_ind= 500150;
# end_ind= 533334;

start_ind_str="{f_ind:010d}".format(f_ind=int(start_ind))
end_ind_str="{f_ind:010d}".format(f_ind=int(end_ind))

print("Case name "+casename+" is read.")
ext_savepath=home_dir+'/postproc/results/ExtSurf/'

# path_UV=home_dir+'/data/MITgcm/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
path_UV=home_dir+'/postproc/results/ExtSurf/'+groupname+'/'+casename+'/'
print("Full path of case extraction: "+path_UV)

path_rawfiles=home_dir+'/scratch/MITgcm_cases/'+groupname+'/'+casename+'/'
print("Raw file is stored at: " + path_rawfiles)

save_img_path=home_dir+'/postproc/img/'
if not os.path.exists(save_img_path):
        os.makedirs(save_img_path)

save_img_path=home_dir+'/postproc/img/'
save_tsr_path=home_dir+'/postproc/tsr/'

#%%
dx=500
dy=500
Lx=480e3
Ly=960e3
xc=np.linspace(dx/2,Lx-dx/2,int(Lx//dx))
yc=np.linspace(dy/2,Ly-dy/2,int(Ly//dy))


# sampling frequency
dt_smp=0.125;
rho_0=1034
freq_f=1.27353
#%%
zeta_surf=Exf_UVP_2D(dx,dy,Lx,Ly,'U'+loc_str,'V'+loc_str,'P'+loc_str,ext_savepath,groupname,
                     casename,start_ind,end_ind)
zeta_surf.load_U()
zeta_surf.load_V()
zeta_surf.calc_rel_vorticity()
zeta_surf.calc_rotary_spectrum() 

U_start_ind_str=zeta_surf.U_fileind[0]
U_end_ind_str=zeta_surf.U_fileind[-1]
V_start_ind_str=zeta_surf.V_fileind[0]
V_end_ind_str=zeta_surf.V_fileind[-1]
#%%
plt.figure(figsize=(3, 4))
plt.pcolormesh(xc/1e3,yc/1e3,zeta_surf.vort[-1]/zeta_surf.f_2D,cmap='bwr');
plt.axis('scaled')
plt.xlabel('x [km]');plt.ylabel('y [km]')
plt.title(r'Relative Vorticity $\zeta/f$')
plt.colorbar();plt.clim([-1,1])
plt.savefig(save_img_path+'vort_'+loc_str+'_'+casename+'_'+U_end_ind_str+'.png',dpi=300)
plt.close()
tch.save(zeta_surf.vort,save_tsr_path+'Vort_'+loc_str+'_'+casename+'_'+
         U_start_ind_str+'_'+U_end_ind_str+'.pt',pickle_protocol=4)
delattr(zeta_surf, 'vort')
print('zeta_surf.vort is deleted.')

# #%%
# exf=utils.rdmds('/aos/home/tpeng/postproc/results/ExtSurf/rel_ucon_12/EXF/EXFtau.0000530014')

# plt.pcolormesh(xc/1e3,yc/1e3,exf[1]);plt.axis('scaled');plt.colorbar();plt.title(r'$U_{10} [m/s]$')
# plt.show()
# plt.pcolormesh(xc/1e3,yc/1e3,exf[0]);plt.axis('scaled');plt.colorbar();plt.title(r'$\tau_x [N/m^2]$')
# plt.show()
#%%
def plot_rotary_spectra(zeta,loc_str):
    rotary_spec=tch.mean(zeta.PRR[:,3:-3,3:-3],dim=[1,2]).real
    Nfile=len(zeta.U)
    nt=tch.linspace(0,1/0.125/2,Nfile//2)
    plt.loglog(nt,(rotary_spec[1:Nfile//2+1])),plt.loglog(nt,tch.flipud(rotary_spec[Nfile//2:]),'--'),plt.legend(('positive - CCW','negative - CW'));
    plt.xlabel('cpd'),plt.ylabel(r'$P_{RR}$'),plt.title('Rotary spectrum -- '+loc_str)
    plt.ylim([1e-3,1e3])
    
def plot_rotary_spectra_geo(zeta,loc_str):
    rotary_spec=tch.mean(zeta.PRR_geo[:,3:-3,3:-3],dim=[1,2]).real
    Nfile=len(zeta.ugeo)
    nt=tch.linspace(0,1/0.125/2,Nfile//2)
    plt.loglog(nt,(rotary_spec[1:Nfile//2+1])),plt.loglog(nt,tch.flipud(rotary_spec[Nfile//2:]),'--'),
    plt.legend(('positive - CCW','negative - CW','positive - CCW -- geo','negative - CW -- geo'));
    plt.xlabel('cpd'),plt.ylabel(r'$P_{RR}$'),plt.title('Rotary spectrum -- '+loc_str)
    plt.ylim([1e-3,1e3])


zeta_surf.load_P()
P_start_ind_str=zeta_surf.P_fileind[0]
P_end_ind_str=zeta_surf.P_fileind[-1]

zeta_surf.calc_rotary_spectrum_geo()

plt.figure()
plot_rotary_spectra(zeta_surf,loc_str)
plot_rotary_spectra_geo(zeta_surf,loc_str)
plt.savefig(save_img_path+'rotary_spec'+loc_str+'_'+casename+'_'+
         U_start_ind_str+'_'+U_end_ind_str+'.png',dpi=300)

tch.save(zeta_surf.PRR,save_tsr_path+'PRR_'+loc_str+'_'+casename+'_'+
         U_start_ind_str+'_'+U_end_ind_str+'.pt',pickle_protocol=4)
tch.save(zeta_surf.PRR,save_tsr_path+'PRR_geo_'+loc_str+'_'+casename+'-'+
         P_start_ind_str+'_'+P_end_ind_str+'.pt',pickle_protocol=4)
delattr(zeta_surf, 'PRR'); 
print('zeta_surf.PRR is deleted.')
delattr(zeta_surf, 'PRR_geo'); 
print('zeta_surf.PRR_geo is deleted.')
#%%
niter_max=45;
zeta_surf.calc_vort_eta(niter_max)
tch.save(zeta_surf.Vort_eta,save_tsr_path+'Vort_eta_'+loc_str+'_'+casename+'_'+
         P_start_ind_str+'_'+P_end_ind_str+'.pt',pickle_protocol=4)
delattr(zeta_surf, 'Vort_eta')
print('zeta_surf.Vort_eta is deleted.')

#%%  print results
