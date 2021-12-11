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

os.mkdir(pythonpath + '/Ez3d_' + tstamp)

c = scipy.constants.c
gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
ts = LpaDiagnostics(path2hdf5)

N_iters = len(ts.iterations)
it = ts.iterations[N_iters-1]
#ax1 = plt.subplot(gs[0])
#ax2 = plt.subplot(gs[0])
steps = np.arange(0,220000,20000)
for it in  steps:
        Ex, info_Ex = ts.get_field( iteration=it, field='E', coord='x',
                                        slicing_dir = 'x', slicing=0.0)
        Ey, info_Ey = ts.get_field( iteration=it, field='E', coord='y',
                                        slicing_dir = 'x', slicing=0.0)
        Ez, info_Ez = ts.get_field( iteration=it, field='E', coord='z',
                                        slicing_dir = 'x', slicing=0.0)


        Bx, info_Bx = ts.get_field( iteration=it, field='B', coord='x',
                                        slicing_dir = 'x', slicing=0.0)
        By, info_By = ts.get_field( iteration=it, field='B', coord='y',
                                        slicing_dir = 'y', slicing=0.0)
        #Fx = Ex - c*By
        #Fy = Ey + c*Bx
        shE = np.shape(Ex)
        shB = np.shape(Bx)
        absmax=max(np.max(np.max(np.abs(Ez)))/10,np.abs(np.min(np.min(np.abs(Ez)))))/10
        print(absmax,'E shape: ', shE, '   B shape: ', shB)
        plt.imshow(Ez[:,int(shE[1]/2),:])#,cmap='RdBu',vmin=-1000,vmax=1000)#extent=[info_Ez.zmin*1e3,info_Ez.zmax*1e3,info_Ez.ymin*1e3,info_Ez.ymax*1e3],cmap='RdBu',vmin=-absmax,vmax=absmax,aspect='auto')
        plt.ylabel('$z$ (mm)')
        plt.xlabel('$x$ (mm)')
#        plt.title('Ez')
        plt.savefig(path2hdf5 + 'Ez_3d'+str(it))
        shutil.move(path2hdf5 + 'Ez_3d'+str(it)+'.png',pythonpath + '/Ez3d_' + tstamp) 

