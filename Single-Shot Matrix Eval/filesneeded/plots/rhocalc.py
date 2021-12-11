import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import gridspec
import scipy.constants
import os
import time
import shutil
from openpmd_viewer import OpenPMDTimeSeries
from openpmd_viewer.addons import LpaDiagnostics



#Path from this python script to the HDF5 files
path2hdf5 = '/data1/aard/cphillips/cslabcode/diags/slabmovingwarpx/'

pythonpath = '/data1/aard/cphillips/cslabcode'

tstamp = time.strftime('%Y%m%d_%H%M%S')

os.mkdir(pythonpath + '/rho_' + tstamp)

c = scipy.constants.c
gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
ts = LpaDiagnostics(path2hdf5)

N_iters = len(ts.iterations)
it = ts.iterations[N_iters-1]
#ax1 = plt.subplot(gs[0])
#ax2 = plt.subplot(gs[0])
steps = np.arange(0,220000,20000)
for it in  steps:
        rho, info_rho = ts.get_field( iteration=it, field='rho')
        shrho = np.shape(rho)
        absmax=max(np.max(np.max(np.abs(rho)))/10,np.abs(np.min(np.min(np.abs(rho)))))/10
        print(absmax)

#ax1.imshow(Ez[int(shE[0]/2),:,:],extent=[info_Ez.zmin*1e3,info_Ez.zmax*1e3,info_Ez.ymin*1e3,info_Ez.ymax*1e3],cmap='RdBu',vmin=-absmax,vmax=absmax,aspect='auto')
        plt.imshow(rho[:,int(shrho[1]/2),:])
#Ez[:,int(shE[1]/2),:])#,cmap='RdBu',vmin=-1000,vmax=1000)#extent=[info_Ez.zmin*1e3,info_Ez.zmax*1e3,info_Ez.ymin*1e3,info_Ez.ymax*1e3],cmap='RdBu',vmin=-absmax,vmax=absmax,aspect='auto')
        #plt.ylabel('$z$ (mm)')
        #plt.xlabel('$x$ (mm)')
#        plt.title('Ez')
        plt.savefig(path2hdf5 + 'rho'+str(it))
        shutil.move(path2hdf5 + 'rho'+str(it)+'.png',pythonpath + '/rho_' + tstamp) 

