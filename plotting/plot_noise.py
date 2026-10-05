#Script for generating Fig. 6 of "Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

L=3
lag = 1
data = pd.read_pickle('ERRORS_L='+str(L)+'_lag='+str(lag)+'.pkl')
lag1_angle = data['angle']
lag1_r = data['anisotropy']

lag = 3
data = pd.read_pickle('ERRORS_L='+str(L)+'_lag='+str(lag)+'.pkl')
lag3_angle = data['angle']
lag3_r = data['anisotropy']

lag = 5
data = pd.read_pickle('ERRORS_L='+str(L)+'_lag='+str(lag)+'.pkl')
lag5_angle = data['angle']
lag5_r = data['anisotropy']

noise = [0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
noise = [0,5e-6,5e-5,5e-4,5e-3,5e-2]
cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

errs_means_1 = lag1_r.mean(axis=2)
errs_means_3 = lag3_r.mean(axis=2)
errs_means_5 = lag5_r.mean(axis=2)

fig, axes = plt.subplots(3, 1,figsize=(6.4, 4.8))
size = fig.get_size_inches()
axes[0].text(-0.15, 1.05, 'd)', transform=axes[0].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
axes[1].text(-0.15, 1.05, 'e)', transform=axes[1].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
axes[2].text(-0.15, 1.05, 'f)', transform=axes[2].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
for i in range(4):
    axes[0].plot(noise[1:],errs_means_1[1:,i],'o-',color=cols[i])
    axes[1].plot(noise[1:],errs_means_3[1:,i],'o-',color=cols[i])
    axes[2].plot(noise[1:],errs_means_5[1:,i],'o-',color=cols[i])
axes[0].grid()
axes[0].set_xscale('log')
axes[0].set_ylabel(r'$\langle RMSE(r)\rangle$',size=12)
axes[0].tick_params(axis='y', which='major', labelsize=12)
axes[0].tick_params(axis='x',label1On=False)
#axes[0].set_yticks([0,1,2],['0','1.0','2.0'])
axes[0].set_ylim(0,2.5)
axes[0].text(1e-5,2.1,r'$L=3, \delta=1$', size=12)

axes[1].grid()
axes[1].set_xscale('log')
axes[1].set_ylabel(r'$\langle RMSE(r)\rangle$',size=12)
axes[1].tick_params(axis='y', which='major', labelsize=12)
axes[1].tick_params(axis='x',label1On=False)
#axes[1].set_yticks([0,1,2],['0','1.0','2.0'])
axes[1].set_ylim(0,2.5)
axes[1].text(1e-5,2.1,r'$L=3, \delta=3$', size=12)

axes[2].grid()
axes[2].set_xscale('log')
axes[2].set_ylabel(r'$\langle RMSE(r)\rangle$',size=12)
axes[2].tick_params(axis='y', which='major', labelsize=12)
axes[2].set_xlabel(r'Noise level, $\eta$',size=12)
axes[2].set_ylim(0,2.5)
axes[2].text(1e-5,2.1,r'$L=3, \delta=5$', size=12)
plt.savefig('noise_r.pdf',format='pdf', dpi=600, bbox_inches='tight')
  
plt.show()

errs_means_1 = lag1_angle.mean(axis=2)
errs_means_3 = lag3_angle.mean(axis=2)
errs_means_5 = lag5_angle.mean(axis=2)

fig, axes = plt.subplots(3, 1,figsize=(6.4, 4.8))
axes[0].text(-0.15, 1.05, 'a)', transform=axes[0].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
axes[1].text(-0.15, 1.05, 'b)', transform=axes[1].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
axes[2].text(-0.15, 1.05, 'c)', transform=axes[2].transAxes,
            fontsize=12, fontweight='normal', va='top', ha='left')
size = fig.get_size_inches()
for i in range(4):
    if i  == 0 or i == 3: axes[0].plot(noise[1:-1],errs_means_1[1:-1,i],'o-',color=cols[i])
    else: axes[0].plot(noise[1:],errs_means_1[1:,i],'o-',color=cols[i])
    axes[1].plot(noise[1:],errs_means_3[1:,i],'o-',color=cols[i])
    axes[2].plot(noise[1:],errs_means_5[1:,i],'o-',color=cols[i])
axes[0].grid()
axes[0].set_xscale('log')
axes[0].set_ylabel(r'$\langle RMSE(\theta)\rangle$',size=12)
axes[0].tick_params(axis='y', which='major', labelsize=12)
axes[0].tick_params(axis='x',label1On=False)
#axes[0].set_yticks([0,10,20],)
axes[0].text(1e-5,10,r'$L=3, \delta=1$', size=12)


axes[1].grid()
axes[1].set_xscale('log')
axes[1].set_ylabel(r'$\langle RMSE(\theta)\rangle$',size=12)
axes[1].tick_params(axis='y', which='major', labelsize=12)
axes[1].tick_params(axis='x',label1On=False)
#axes[1].set_yticks([0,1,2],['0','1.0','2.0'])
axes[1].set_ylim(0,14)
axes[1].text(1e-5,10,r'$L=3, \delta=3$', size=12)

axes[2].grid()
axes[2].set_xscale('log')
axes[2].set_ylabel(r'$\langle RMSE(\theta)\rangle$',size=12)
axes[2].tick_params(axis='y', which='major', labelsize=12)
#axes[2].set_xlabel('Noise level',size=12)
axes[2].set_ylim(0,14)
axes[2].tick_params(axis='x',label1On=False)
axes[2].text(1e-5,10,r'$L=3, \delta=5$', size=12)

axes[0].set_ylim(0,14)
axes[0].set_xlim(axes[2].get_xlim())

fig.align_ylabels()

plt.savefig('noise_angle.pdf',format='pdf', dpi=600, bbox_inches='tight')
    
plt.show()