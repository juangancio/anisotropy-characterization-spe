import numpy as np
import pandas as pd
import xarray as xr
import tqdm, os, pickle, math
import multiprocess as mp
import matplotlib.pyplot as plt

#This scrips creates the panels used for Fig.2 of the article
#"Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.


models = ['AP','AP','GLE','CBM']
mods = ['chaos','ulrich','chaos','chaos']
texts = ['a)','b)','c)','d)']
titles = ['Aliev-Panfilov (fibrillation)',
          'Aliev-Panfilov (artificial)',
          'Ginzburg-Landau equation',
          'Cubic Barkley model']

labels = [r'$u(x,y)$',r'$u(x,y)$',r'$\phi(x,y)$',r'$u(x,y)$']

mult = 1

phi = 45
D = 200
k = 0
for model in models:
    mod = mods[k]
    txt = texts[k]
    title = titles[k]
    if k == 2: map = 'hsv'; interp = 'nearest'
    else: map = 'jet'; interp = 'gaussian'
    k += 1

    if model == 'AP':
        if mod == 'chaos':
            D1 = np.array([.3,.4,0.3,.4,.4,.6])*mult
            D2 = np.array([.2,.2,.1,.1,0.08,.1])*mult
        if mod == 'ulrich':
            D1 = np.array([ 1.5,2, 1.8, 2,2,2.4])*mult
            D2 = np.array([1,1, .6, .5, .4, .4])*mult
        if mod == 'spiral':
            D1 = np.array([ 1.4, 1.8, 2,2])*mult
            D2 = np.array([0.7, .6, .5, .4])*mult

    if model == 'CBM':
        D1 = np.array([6,8,9,8,10,12])*mult
        D2 = np.array([4,4,3,2,2,2])*mult
    if model == 'GLE':
        D1 =  np.array([.6,.8,.9,.8,1,1.2])*mult
        D2 =  np.array([.4,.4,.3,.2,.2,.2])*mult


    d1 = D1[4]; d2 = D2[4]
    if np.mod(d1,1) == 0: d1 = int(d1)
    if np.mod(d2,1) == 0: d2 = int(d2)

    if os.path.isfile('/Volumes/T7/RD_data/'+model+'/'+model+'_d1_'
            +str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi))+'_'+mod+'.h5'):
        file = xr.load_dataset('/Volumes/T7/RD_data/'+model+'/'+model+'_d1_'
                +str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi))+'_'+mod+'.h5', engine="h5netcdf")
        data = file["data"].to_numpy()
    elif os.path.isfile('/Volumes/T7/RD_data/'+model+'/'+model+'_d1_'
            +str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi))+'_'+mod+'.csv'):
        data = np.array(
            pd.read_csv('/Volumes/T7/RD_data/'+model+'/'+model+'_d1_'
            +str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi))+'_'+mod+'.csv')
            )
    else: raise Exception("No valid file")

    fig, axes = plt.subplots(1, 1,figsize=(6.4, 5.5))
    plt.imshow(data[0].reshape((D,D)),origin='lower',cmap=map,aspect='auto',interpolation=interp)
    plt.xlabel(r'$x$',size=20); plt.ylabel(r'$y$',size=20,rotation=90)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=20)
    cb.set_label(labels[k-1], rotation=270,fontsize=20, labelpad=25)
    plt.tick_params(axis='both',label1On=False)
    plt.title(title,size=20)
    axes.text(-0.08, 1.07, txt, transform=axes.transAxes,
            fontsize=20, fontweight='normal', va='top', ha='left')
    plt.savefig('snapshot_'+txt[0]+'.pdf',format='pdf', dpi=600, bbox_inches='tight')
    plt.show()