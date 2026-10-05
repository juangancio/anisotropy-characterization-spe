## Script for the spatial analysis of RD systems.
# It calculates the Spatial Permutation Entropy (SPE) for each time and in both axes direction,
# and outputs their mean and standard deviation. The script is prepared to load files with 
# different orientation of the diffusion axes, such as they swept a certain angle range.
# The script performs the analysis, and plots and saves the results.
#
# The script now supports the calculation of "amplitude aware" SPE, this is controlled by the
# 'mode' parameter (mode='usual' for usual SPE, or mode='amplitude' for AASPE).

import numpy as np
import pandas as pd
import xarray as xr
import tqdm, os, pickle, math
import multiprocess as mp
import matplotlib.pyplot as plt

from spe_utils import (
    h_entropy,
    v_entropy,
    perm_indices,
    entropy,
    diag1_entropy,
    diag2_entropy
)

def autocorr(x,lags):

    mean=np.mean(x)
    var=np.var(x)
    xp=x-mean
    corr=np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

    return corr[:lags]

def h_auto(data):
    code = []
    for i in range(len(data[:,0])):
        code.append(autocorr(data[i,:],2)[1])
    return np.mean(code)

def v_auto(data):
    code = []
    for i in range(len(data[0,:])):
        code.append(autocorr(data[:,i],2)[1])
    return np.mean(code)

def weighting(ts,wl,lag,A):
    w = []
    for i in range(len(ts) - (wl - 1)*lag):
        word = []
        for j in range(0,(wl )*lag,lag):
            word.append(ts[i+j])
        w.append(A * np.mean(np.abs(word)) + (1 - A) * np.mean(np.abs(np.diff(word))))
    return w

def w_probabilities(code,weights,L):
    get_indexes = lambda x, xs: [k for (y, k) in zip(xs, range(len(xs))) if x == y]
    probs =[]
    norm = sum(weights)
    for i in range(1,math.factorial(L)+1):
        
        probs=probs + [sum(weights[get_indexes(i,code)])/norm
                   ]
        
    return probs

def h_aa_entropy(data,L,lag):
    code=[]
    weights=[]
    for i in range(len(data[:,0])):
        code.append(perm_indices(data[i,:],L,lag))
        weights.append(weighting(data[i,:],L,lag,.5))
    return entropy(w_probabilities(np.array(code).flatten(),np.array(weights).flatten(),L))/np.log(math.factorial(L))

def v_aa_entropy(data,L,lag):
    code=[]
    weights=[]
    for i in range(len(data[0,:])):
        code.append(perm_indices(data[:,i],L,lag))
        weights.append(weighting(data[:,i],L,lag,.5))
    return entropy(w_probabilities(np.array(code).flatten(),np.array(weights).flatten(),L))/np.log(math.factorial(L))

def run_h(t):
    return h_entropy(data[t,:].reshape((D,D)),L,lag)
def run_v(t):
    return v_entropy(data[t,:].reshape((D,D)),L,lag)

def run_diag1(t):
    return diag1_entropy(data[t,:].reshape((D,D)),L,lag)
def run_diag2(t):
    return diag2_entropy(data[t,:].reshape((D,D)),L,lag)

def run_aa_h(t):
    return h_aa_entropy(data[t,:].reshape((D,D)),L,lag)
def run_aa_v(t):
    return v_aa_entropy(data[t,:].reshape((D,D)),L,lag)

def run_auto_h(t):
    return h_auto(data[t,:].reshape((D,D)))
def run_auto_v(t):
    return v_auto(data[t,:].reshape((D,D)))

#################################################################################

L = 3 #length of the symbols
lag = 1 #spatial lag of the symbols



angles = np.linspace(0,90,19) #range of angles to be scanned
cores = 8 #number of cores to use in the computation

mode = 'usual' #'usual', 'diagonal' or 'amplitude'
#model = 'GLE'  #'CBM', 'AP', or 'GLE'
#mod = 'chaos'
#################################################################################


#AP chaos
'''D1 = [.3,.4,0.3,.4,.4,.6]
D2 = [.2,.2,.1,.1,0.08,.1]

D1 = np.array(D1)/2
D2 = np.array(D2)/2


D1 = [.3,.6]
D2 = [.2,.1]
'''
#D1 = [.4,0.3,.4,.4]
#D2 = [.2,.1,.1,0.08]

#AP ulrich
'''D1 = [ 2, 1.8, 2,2]
D2 = [1, .6, .5, .4]

D1 = [ 1.5,2.4]
D2 = [1, .4]

D1 = [ 1.5,2, 1.8, 2,2,2.4]
D2 = [1,1, .6, .5, .4, .4]

D1 = [ 0.75,1, 0.9, 1,1,1.2]
D2 = [0.5,0.5, 0.3, 0.25, 0.2, 0.2]'''

#GLE
'''D1 = [.8,.9,.8,1]
D2 = [.4,.3,.2,.2]

D1 = [.6,1.2]
D2 = [.4,.2]

D1 = [.6,.8,.9,.8,1,1.2]
D2 = [.4,.4,.3,.2,.2,.2]

D1 = np.array(D1)*2
D2 = np.array(D2)*2
'''

#CBM

'''D1 = [8,9,8,10]
D2 = [4,3,2,2]

D1 = [6,12]
D2 = [4,2]

D1 = [6,8,9,8,10,12]
D2 = [4,4,3,2,2,2]
D1 = np.array(D1)*2
D2 = np.array(D2)*2'''

mods = ['chaos','ulrich','chaos','chaos']
for noise_level in [0]:#,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]:
    
    #noise_level = 0.000001
    for mult in [1]:
        k = 0
        for model in ['AP']:#,'AP','CBM','GLE']:  
            mod = mods[k]
            k += 1
            if mult == 2 and mod == 'ulrich': continue
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


            for i in range(len(D1)):

                d1 = D1[i] # diffusion coefficients
                d2 = D2[i] # in each direction

                if np.mod(d1,1) == 0: d1 = int(d1)
                if np.mod(d2,1) == 0: d2 = int(d2)

                H=[]; s_H=[]; V=[]; s_V=[]
                Ha=[]; s_Ha=[]; Va=[]; s_Va=[]

                for phi in tqdm.tqdm([60]):#angles):

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

                    #data = data + np.random.normal(scale=np.max(abs(data))*noise_level,size=data.shape)

                    h = []; v = []
                    D = int(np.sqrt(data.shape[1]))

                    if mode == 'usual':
                        if __name__ == '__main__':
                            pool = mp.Pool(cores)       
                            #h_code=pool.map(get_code,range(data.shape[0]))
                            h = pool.map(run_h,range(data.shape[0]))
                            v = pool.map(run_v,range(data.shape[0]))
                            ha = pool.map(run_auto_h,range(data.shape[0]))
                            va = pool.map(run_auto_v,range(data.shape[0]))
                            pool.close()  
                            pool.join() 

                    elif mode == 'diagonal':
                        if __name__ == '__main__':
                            pool = mp.Pool(cores)       
                            #h_code=pool.map(get_code,range(data.shape[0]))
                            h = pool.map(run_diag1,range(data.shape[0]))
                            v = pool.map(run_diag2,range(data.shape[0]))
                            pool.close()  
                            pool.join() 

                    elif mode == 'amplitude':
                        if __name__ == '__main__':
                            pool = mp.Pool(cores)       
                            #h_code=pool.map(get_code,range(data.shape[0]))
                            h = pool.map(run_aa_h,range(data.shape[0]))
                            v = pool.map(run_aa_v,range(data.shape[0]))
                            pool.close()  
                            pool.join() 
                    else:
                        raise Exception("Specified mode is not valid, try 'usual' or 'amplitude'.")

                    H.append(np.mean(h)); V.append(np.mean(v))
                    s_H.append(np.std(h)); s_V.append(np.std(v))

                    #Ha.append(np.mean(ha)); Va.append(np.mean(va))
                    #s_Ha.append(np.std(ha)); s_Va.append(np.std(va))

                H = np.array(H); V = np.array(V); s_H = np.array(s_H); s_V = np.array(s_V)

                #Ha = np.array(Ha); Va = np.array(Va); s_Ha = np.array(s_Ha); s_Va = np.array(s_Va)


                #est = H / V
                #err = np.sqrt((s_H/V)**2 + (H/(V**2) * s_V)**2)
                #fit = np.polyfit(angles, est, deg=1)
                #fit3 = np.polyfit(angles, est, deg=3)


                #est_a = Va / Ha
                #err_A = np.sqrt((s_Va/Ha)**2 + (Va/(Ha**2) * s_Va)**2)

                #data_to_save = { "H_h":H, "H_v":V, "H_h_std":s_H, "H_v_std":s_V, "ACF_h":Ha, "ACF:_v":Va, "sACF_h":s_Ha, "sACF:_v":s_Va}
                
                data_to_save = { "H_h":H, "H_v":V, "H_h_std":s_H, "H_v_std":s_V,}
                
                '''if mode == 'usual':
                    with open('processed/spe_axes_d1_'+str(d1)+'_d2_'+str(d2)+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl', "wb") as f:
                            pickle.dump(data_to_save, f)
                if mode == 'diagonal':
                    with open('processed/spe_axes_diag_d1_'+str(d1)+'_d2_'+str(d2)+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl', "wb") as f:
                            pickle.dump(data_to_save, f)
                elif mode == 'amplitude':
                    with open('angle_test/aa_spe_axes_d1_'+str(d1)+'_d2_'+str(d2)+'_L_'+str(L)+'_lag_'+str(lag)+'_chaos.pkl', "wb") as f:
                            pickle.dump(data_to_save, f)'''




'''fig, axes = plt.subplots(1, 1,figsize=(6.4, 5.5))
plt.imshow(data[0].reshape((D,D)),origin='lower',cmap='plasma',aspect='auto')
plt.tick_params(axis='both',label1On=False)
axes.text(-0.15, 1.05, r'a)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.savefig('spiral_snapshot.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()'''


'''fig, axes = plt.subplots(1, 1,figsize=(6.4, 5.5))
toplot=data[50].reshape((D,D))
plt.contour(toplot[40:,:-40],origin='lower',levels=[.5],aspect='square',colors='k',linestyles='dashed')
plt.tick_params(axis='both',label1On=False)
axes.text(-0.15, 1.05, r'a)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.savefig('spiral_snapshot_methods.png',format='png', dpi=600, bbox_inches='tight')
plt.show()'''