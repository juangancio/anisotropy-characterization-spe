#Script for generating Fig. 8 of "Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as op 

B = []; 
cols = ['#1f77b4',
 '#ff7f0e',
 '#2ca02c',
 '#d62728',
 '#9467bd',
 '#8c564b']
fig, axes = plt.subplots(1, 1)
for L in [3,4,5]:
    raw_data = pd.read_csv('exponents_L='+str(L)+'_lag=1.csv')
    i = 0
    B.append(raw_data['b'].mean())
    '''for r in [1.5,2,3,4,5,6]:
        filt_data = raw_data[raw_data['r'] == r]'''
    for r in ['AP_chaos','AP_ulrich','GLE','CBM']:
        filt_data = raw_data[raw_data['model'] == r]
        b = -filt_data['b']
        p = plt.plot([L]*len(b),b,'o',color = cols[i],alpha=.3)
        
        i += 1
    i = 0
    for r in ['AP_chaos','AP_ulrich','GLE','CBM']:
        filt_data = raw_data[raw_data['model'] == r]
        b = -filt_data['b']
        #p = plt.plot([L]*len(b),b,'o',color = cols[i],alpha=.3)
        plt.plot(L+.07,np.mean(b),'<',color = cols[i],markersize=8)
        i += 1
        

plt.xticks([3,4,5])
plt.text(4.6,-0.3,r'$\mathbf{(\delta=1)}$',size=14)
plt.xlabel(r'Word length, $L$',size=14)
plt.ylabel(r'$b$',size=14)
plt.tick_params(axis='both', which='major', labelsize=14)
axes.text(-0.2, 1.05, 'a)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.grid()
plt.savefig('L_effect.pdf',format='pdf', bbox_inches='tight')
plt.show()
#filt_data = raw_data[raw_data['model'] != 'AP_chaos']
#filt_data = raw_data
#filt_data = filt_data[filt_data['d1']/filt_data['d2']  == 2]
#x, y = filt_data['angle'], filt_data['H']/filt_data['V']
fig, axes = plt.subplots(1, 1)
L = 3
labels = ['AP (fib.)', 'AP (art.)', 'GL', 'CB']
for lag in [1,2,3,4,5]:
    raw_data = pd.read_csv('exponents_L='+str(L)+'_lag='+str(lag)+'.csv')
    i = 0
    '''for r in [1.5,2,3,4,5,6]:
        filt_data = raw_data[raw_data['r'] == r]'''
    for r in ['AP_chaos','AP_ulrich','GLE','CBM']:
        filt_data = raw_data[raw_data['model'] == r]
        b = -filt_data['b']
        p = plt.plot([lag]*len(b),b,'o',color = cols[i],alpha=.3)
        #if lag == 1:
        #    p = plt.plot([lag]*len(b),b,'o',color = cols[i],alpha=.25,label=labels[i])
        i += 1
    i = 0
    for r in ['AP_chaos','AP_ulrich','GLE','CBM']:
        filt_data = raw_data[raw_data['model'] == r]
        b = -filt_data['b']
        #p = plt.plot([L]*len(b),b,'o',color = cols[i],alpha=.3)
        plt.plot(lag+.15,np.mean(b),'<',color = cols[i],markersize=8)
        i += 1
plt.xticks([1,2,3,4,5])
plt.xlabel(r'Spatial lag, $\delta$',size=14)
plt.ylabel(r'$b$',size=14)
plt.text(4.2,-.175,r'$\mathbf{(L=3)}$',size=14)
plt.tick_params(axis='both', which='major', labelsize=14)
axes.text(-0.2, 1.05, 'b)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.grid()
plt.savefig('lag_effect.pdf',format='pdf', bbox_inches='tight')

'''legend = axes.legend(fontsize=14,loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=4)
def export_legend(legend, filename="legend_2.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)
export_legend(legend)'''

plt.show()