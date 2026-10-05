import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import scipy.optimize as op 

## Script for plotting the analized data obtained from the different system.
# It plots the majority of the figures included in supplementary material,
# and Fig. 7b of the main text.
# 
# Figures included in "Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.
 

def curve(angle,r):
    angle = np.radians(angle)
    return (np.cos(angle)**2+r*np.sin(angle)**2+(r-1)*np.cos(angle)*np.sin(angle))/(r*np.cos(angle)**2+np.sin(angle)**2+(r-1)*np.cos(angle)*np.sin(angle))
    #return (np.cos(angle)**2+r*np.sin(angle)**2)/(r*np.cos(angle)**2+np.sin(angle)**2)
    #return (np.cos(angle)**2+r*np.sin(angle)**2+2*(r-1)*np.cos(angle)*np.sin(angle))/(r*np.cos(angle)**2+np.sin(angle)**2+2*(r-1)*np.cos(angle)*np.sin(angle))
    #return np.sqrt(np.cos(angle)**2+r*np.sin(angle)**2)/np.sqrt(r*np.cos(angle)**2+np.sin(angle)**2)
def curve_2(angle,a,b,c):
    angle = np.radians(angle)
    return (b*np.cos(angle)**2+a*np.sin(angle)**2+c*np.cos(angle)*np.sin(angle))/(a*np.cos(angle)**2+b*np.sin(angle)**2+c*np.cos(angle)*np.sin(angle))
    #return (np.cos(angle)**2+r*np.sin(angle)**2+2*(r-1)*np.cos(angle)*np.sin(angle))/(r*np.cos(angle)**2+np.sin(angle)**2+2*(r-1)*np.cos(angle)*np.sin(angle))
    #return np.sqrt(np.cos(angle)**2+r*np.sin(angle)**2)/np.sqrt(r*np.cos(angle)**2+np.sin(angle)**2)
def curve_14(angle,r):
    angle = np.radians(angle)
    return (np.cos(angle)**2+r*np.sin(angle)**2)/(r*np.cos(angle)**2+np.sin(angle)**2)

def curve_H(angle,a,b):
    angle = np.radians(angle)
    return 1/(a*np.cos(angle)**2+b*np.sin(angle)**2)

def curve_H_2(angle,a,b,c):
    angle = np.radians(angle)
    return 1/(a*np.cos(angle)**2+b*np.sin(angle)**2+c*np.cos(angle)*np.sin(angle))

def curve_new(angle,r):
    angle = np.radians(angle)
    return (np.cos(angle)**2+r*np.sin(angle)**2)/(r*np.cos(angle)**2+np.sin(angle)**2)**0.2

def calculate_elements(angle,d1,d2):
    angle = np.radians(angle)
    a = d1*np.cos(angle)**2+d2*np.sin(angle)**2
    b = d2*np.cos(angle)**2+d1*np.sin(angle)**2
    c = (d1-d2)*np.cos(angle)*np.sin(angle)
    return a,b,c

def curve_rate(X,c):
    d11,d22,d12 = X
    return ((d22 / (d11)) ** c)
def curve_hor(X,b,c):
    d11,d22,d12 = X
    return (b*(d11) ** c)
def curve_ver(X,b,c):
    d11,d22,d12 = X
    return (b*(d22) ** c)
def curve_r(X,b):
    angle,r = X
    r = r[0]
    angle = np.radians(angle)
    return ((np.cos(angle)**2+r*np.sin(angle)**2)/(r*np.cos(angle)**2+np.sin(angle)**2))**b

def curve_hor_log(X,b,c):
    d11,d22,d12 = X
    return b*np.log(d11)+c
def curve_ver_log(X,b,c):
    d11,d22,d12 = X
    return b*np.log(d22)+c

cycol = cycle('bgrcmk')
angles=np.linspace(0,90,19)
d1=0.08
d2=0.02


L=3
lag=4

model = 'AP'
mod='chaos'

ani = ['r=2','r=3','r=4','r=5']
d1 = [.4,0.3,.4,.4]
d2 = [.2,.1,.1,0.08]
cols=[]

d1_to_save = []
d2_to_save = []
angles_to_save = []
entropy_to_save = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p=plt.plot(angles, est,label=ani[i],lw=3)
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_rate,(D11,D22,D12),est)
    if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    cols.append(p[0].get_color())
    d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)


plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.55,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_fib_h.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    s_H = data['H_h_std']
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_hor,(D11,D22,D12),H)
    plt.plot(angles, H,label=ani[i])
    plt.fill_between(angles, H-s_H, H+s_H,alpha=0.5)
    plt.plot(angles,curve_hor((D11,D22,D12),fit[0][0],fit[0][1]),'k--')
    
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'horizontal entropy, $H_h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.62,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_fib_H_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')

    V = data['H_v']
    s_V = data['H_v_std']
    plt.plot(angles, V,label=ani[i])
    plt.fill_between(angles, V-s_V, V+s_V,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'vertical entropy, $H_v$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.62,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_fib_H_v.png',format='png', dpi=600, bbox_inches='tight')
plt.show()



for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
    
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-15,1.105,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_fib_alpha.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


model = 'AP'
mod='ulrich'

ani = ['r=2','r=3','r=4','r=5']
d1 = [ 2, 1.8, 2,2]
d2 = [1, .6, .5, .4]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],lw=3)
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_rate,(D11,D22,D12),est)
    if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.45,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_art_h.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    s_H = data['H_h_std']
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_hor,(D11,D22,D12),H)
    plt.plot(angles, H,label=ani[i])
    plt.fill_between(angles, H-s_H, H+s_H,alpha=0.5)
    plt.plot(angles,curve_hor((D11,D22,D12),fit[0][0],fit[0][1]),'k--')
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'horizontal entropy, $H_h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.61,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_art_H_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')

    V = data['H_v']
    s_V = data['H_v_std']
    plt.plot(angles, V,label=ani[i])
    plt.fill_between(angles, V-s_V, V+s_V,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'vertical entropy, $H_v$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.61,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_art_H_v.png',format='png', dpi=600, bbox_inches='tight')
plt.show()


for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.055,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('AP_art_alpha.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


model = 'GLE'
mod='chaos'

ani = ['r=2','r=3','r=4','r=5']
d1 = [.8,.9,.8,1]
d2 = [.4,.3,.2,.2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],lw=3)
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_rate,(D11,D22,D12),est)
    if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.5,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('GLE_h.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    s_H = data['H_h_std']

    plt.plot(angles, H,label=ani[i])
    plt.fill_between(angles, H-s_H, H+s_H,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'horizontal entropy, $H_h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.92,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('GLE_H_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')

    V = data['H_v']
    s_V = data['H_v_std']
    plt.plot(angles, V,label=ani[i])
    plt.fill_between(angles, V-s_V, V+s_V,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'vertical entropy, $H_v$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.92,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('GLE_H_v.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.6,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('GLE_alpha.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


model = 'CBM'
mod='chaos'

ani = ['r=2','r=3','r=4','r=5']
d1 = [8,9,8,10]
d2 = [4,3,2,2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],lw=3)
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_rate,(D11,D22,D12),est)
    if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.45,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('CBM_h.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

data_to_save = np.c_[np.array(d1_to_save), d2_to_save, angles_to_save, entropy_to_save]
data_to_save = np.array(data_to_save)
data_to_save = pd.DataFrame(data_to_save,columns=['d1','d2','angle','H'])
#data_to_save.to_csv('/Users/juan/Downloads/rguimera-machine-scientist-c790f7cf35b5/test.csv')


for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['H_h']
    s_H = data['H_h_std']

    plt.plot(angles, H,label=ani[i])
    plt.fill_between(angles, H-s_H, H+s_H,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'horizontal entropy, $H_h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.54,r'g)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('CBM_H_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')

    V = data['H_v']
    s_V = data['H_v_std']
    plt.plot(angles, V,label=ani[i])
    plt.fill_between(angles, V-s_V, V+s_V,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'vertical entropy, $H_v$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,.54,r'h)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('CBM_H_v.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.04,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('CBM_alpha.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()



models = ['AP','AP','GLE','CBM']
mods=['chaos','ulrich','chaos','chaos']

ani = ['AP (fib.)','AP (artificial)','GL','CB']

fits = []
fits_h = []; fits_v = []
fits_h_2 = []; fits_v_2 = []

fits_rate = []; 
fits_rate_2 = []; 

rss = []; rss_v1 = []; rss_v2 = []; rss_v3 = []; rss_v4 = []
#############################################################
r = 1.5; fit_r = []
d1 = [.3,1.5,.6,6]
d2 = [.2,1,.4,4]
cols = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3)
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits.append(op.curve_fit(curve,angles,est)[0])
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    cols.append(p[0].get_color())# 'b'


angle = np.linspace(0,1,110)*np.pi/2
#plt.plot(np.degrees(angle),(np.cos(angle)**2+r*np.sin(angle)**2+(r-1)*np.cos(angle)*np.sin(angle))/(r*np.cos(angle)**2+np.sin(angle)**2+(r-1)*np.cos(angle)*np.sin(angle)),'k-',linewidth=3,label='Fit (Eq. 20)')
#plt.plot(np.degrees(angle),(np.cos(angle)**2+r*np.sin(angle)**2)/(r*np.cos(angle)**2+np.sin(angle)**2),'--',color='lightgrey',linewidth=3,label='Fit (Eq. 14)')

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.85,1.33))
plt.text(-10,1.32,r'a)',size=16)
plt.legend(prop={'size': 16},loc=2); 
plt.savefig('h_systems_r=1.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()
#############################################################

r =2; fit_r = []
d1 = [.4,2,.8,8]
d2 = [.2,1,.4,4]
cols = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3)
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    cols.append(p[0].get_color())# 'b'

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.8,1.27))
plt.text(-10,1.26,r'b)',size=16)
plt.legend(prop={'size': 16},loc=2); 
plt.savefig('h_systems_r=2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

r = 3; fit_r = []
d1 = [.3,1.8,.9,9]
d2 = [.1,.6,.3,3]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.37,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=3.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

r = 4; fit_r = []
d1 = [.4,2,.8,8]
d2 = [.1,.5,.2,2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.48,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=4.pdf',format='png', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.4,2,1,10]
d2 = [.08,.4,.2,2]
r = 5; fit_r = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fit_h = op.curve_fit(curve_H_2,angles,H)[0]
    fit_v = op.curve_fit(curve_H_2,angles,V)[0]
    fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.5,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

r = 6; fit_r = []
d1 = [.6,2.4,1.2,12]
d2 = [.1,.4,.2,2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.65,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=6.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()



r = 6; fit_r = []
d1 = [.6,2.4,1.2,12]
d2 = [.1,.4,.2,2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    #fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    #fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    #fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    #fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    #fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit',linewidth=2)

r = 1.5; fit_r = []
d1 = [.3,1.5,.6,6]
d2 = [.2,1,.4,4]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est, linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    #fits.append(op.curve_fit(curve,angles,est)[0])
    #fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    #fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    #fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    #fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    #cols.append(p[0].get_color())# 'b'

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',linewidth=2)

plt.text(65,0.95,r'$r=1.5$',size=16)
plt.text(65,1.44,r'$r=6$',size=16)
plt.text(60,.6,r'$\mathbf{L=4, \delta=1}$',size=16)
plt.xlabel(r'Angle $\theta$ ($^o$)',size=16); plt.ylabel(r'Entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
#plt.text(-10,1.65,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=1.5+6.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()




r = 5; fit_r = []
d1 = [.4]#,2,1,10]
d2 = [.08]#,.4,.2,2]

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    selection = [5]
    plt.plot(angles[selection], est[selection],'--',label=r'Mean entropy ratio measurements, $\langle h_{\hat{v}}\rangle_t$, $\langle h_{\hat{u}}\rangle_t$', linewidth=3, color = 'b')
    plt.plot(angles[selection], est[selection],'--',label=r'Inferred angles', linewidth=3, color = cols[3],markerfacecolor=cols[3])

    #if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.5,label='Std. dev. in time')
    #else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    #fits.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    #fits_h.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    #fits_v.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    #fits_rate.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    #fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)+0.2),'-',color='k',label=r'Model, $\tilde{h}$ (Eq. (9))',linewidth=2.5)
#selection = [5,14]
#plt.plot(angles[selection], est[selection],'--',label=r'Inferred ', linewidth=3, color = cols[3],markerfacecolor=cols[3])

'''plt.text(65,0.95,r'$r=1.5$',size=16)
plt.text(65,1.44,r'$r=6$',size=16)'''

plt.xlabel(r'Angle $\theta$ ($^o$)',size=12); plt.ylabel(r'Entropy ratio, $h$',size=12 )
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(); plt.xlim((0,90))
#plt.text(-10,1.65,r'f)',size=16)
plt.legend(prop={'size': 12}); 
plt.savefig('h_systems_r=5_methods.png',format='png', dpi=600, bbox_inches='tight')
plt.show()






#r = 1.5
d1 = [.3,1.5,.6,6]
d2 = [.2,1,.4,4]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.11,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=1.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


d1 = [.4,2,.8,8]
d2 = [.2,1,.4,4]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.17,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.3,1.8,.9,9]
d2 = [.1,.6,.3,3]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',linewidth=3,label=ani[i], color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.32,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=3.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.4,2,.8,8]
d2 = [.1,.5,.2,2]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.5,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=4.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.4,2,1,10]
d2 = [.08,.4,.2,2]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.55,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.6,2.4,1.2,12]
d2 = [.1,.4,.2,2]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.61,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=6.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

N=200

angles = np.linspace(0,90,N)
mesh = np.meshgrid(angles,angles)
angles = np.sqrt(mesh[0]**2+mesh[1]**2)
angles = 90 * angles/np.sqrt(2*90**2)
s_phi =  np.radians(angles)
plt.imshow(angles,origin='lower')
plt.set_cmap('plasma')
plt.xlabel('horizontal coordinate',size=16)
plt.ylabel('vertical coordinate',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
cbar=plt.colorbar()
cbar.set_label(r'angle, $\theta$',size=16)
cbar.ax.tick_params(labelsize=16) 
#plt.savefig('angle_color.png')
plt.show()

plt.streamplot(np.linspace(0,N,N),np.linspace(0,N,N),np.cos(s_phi),np.sin(s_phi))
plt.xlabel('horizontal coordinate',size=16)
plt.ylabel('vertical coordinate',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
#plt.savefig('angle_lines.png')
plt.show()




angles=np.linspace(0,90,19)
models = ['AP','AP','GLE','CBM']
mods=['chaos','ulrich','chaos','chaos']

ani = ['AP (fib.)','AP (artificial)','GL','CB']

fits_2 = []
fits_h_ii = []; fits_v_ii = []
fits_h_2_ii = []; fits_v_2_ii = []
fits_rate_ii = []; 

rss_2 = []; rss_v1_2 = []; rss_v2_2 = []; rss_v3_2 = []; rss_v4_2= []
#############################################################
r = 1.5; fit_r = []

d1 = [.15,.75,.3,3]
d2 = [.1,0.5,.2,2]

cols = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3)
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    cols.append(p[0].get_color())# 'b'

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.85,1.33))
plt.text(-10,1.32,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=1.5_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()
#############################################################


r = 2; fit_r = []
d1 = [.2,1,.4,4]
d2 = [.1,.5,.2,2]
cols = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3)
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    cols.append(p[0].get_color())# 'b'

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.8,1.27))
plt.text(-10,1.26,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=2_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.15,.9,.45,4.5]
d2 = [.05,.3,.15,1.5]
r = 3; fit_r = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.47,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=3_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.2,1,.4,4]
d2 = [0.05,.25,.1,1]
r = 4; fit_r = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.69,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=4_0.5.pdf',format='png', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.2,1,.5,5]
d2 = [.04,.2,.1,1]
r = 5; fit_r = []

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])


plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.69,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=5_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


d1 = [.3,1.2,.6,6]
d2 = [.05,.2,.1,1]
r = 6; fit_r = []
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8)
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5)
    fits_2.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_ii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_ii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_ii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])


plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.67,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=6_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()





#############################################################
r = 1.5; fit_r = []

d1 = [.6,3,1.2,12]
d2 = [.4,2,.8,8]

fits_3 = []
fits_h_iii = []; fits_v_iii = []
fits_h_2_iii = []; fits_v_2_iii = []
fits_rate_iii = []; 

rss_3 = []; rss_v1_3 = []; rss_v2_3 = []; rss_v3_3 = []; rss_v4_3 = []

for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3,color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

    #type(p) # list
    #type(p[0]) # <class 'matplotlib.lines.Line2D'>
    #cols.append(p[0].get_color())# 'b'

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.85,1.33))
plt.text(-10,1.32,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=1.5_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()
#############################################################


r = 2; fit_r = []
d1 = [.8,4,1.6,16]
d2 = [.4,2,.8,8]

for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p = plt.plot(angles, est,label=ani[i], linewidth=3,color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90)); plt.ylim((0.8,1.27))
plt.text(-10,1.26,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=2_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.6,3.6,1.8,18]
d2 = [.2,1.2,.6,6]
r = 3; fit_r = []
for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.37,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=3_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.8,4,1.6,16]
d2 = [0.2,1,.4,4]
r = 4; fit_r = []
for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.42,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=4_2.pdf',format='png', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.8,4,2,20]
d2 = [.16,.8,.4,4]
r = 5; fit_r = []

for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fit_h = op.curve_fit(curve_H_2,angles,H)[0]
    fit_v = op.curve_fit(curve_H_2,angles,V)[0]
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.62,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=5_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


d1 = [1.2,4.8,2.4,24]
d2 = [.2,.8,.4,4]
r = 6; fit_r = []

for i in range(len(d1)):
    if i ==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i], linewidth=3, color = cols[i])
    if i == 0: plt.fill_between(angles, est-err, est+err,alpha=0.8, color = cols[i])
    else: plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
    fits_3.append(op.curve_fit(curve,angles,est)[0])
    D11,D22,D12 = calculate_elements(angles,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fits_h_iii.append(op.curve_fit(curve_hor,(D11,D22,D12),H)[0])
    fits_v_iii.append(op.curve_fit(curve_ver,(D11,D22,D12),V)[0])
    fits_rate_iii.append(op.curve_fit(curve_rate,(D11,D22,D12),H/V)[0])
    fit_r.append(op.curve_fit(curve_r,(angles,r*np.ones(angles.shape)),est)[0])

plt.plot(angles,curve_r((angles,r*np.ones(angles.shape)),np.mean(fit_r)),'--',color='k',label=r'Fit ($|b|='+str(np.round(np.mean(fit_r),3))+r'$)',linewidth=2)

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.67,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('h_systems_r=6_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()



#r = 1.5
d1 = [.15,0.75,.3,3]
d2 = [.1,.5,.2,2]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.22,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=1.5_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


d1 = [.2,1,.4,4]
d2 = [.1,.5,.2,2]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.33,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=2_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.15,.9,.45,4.5]
d2 = [.05,.3,.15,1.5]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',linewidth=3,label=ani[i], color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.66,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=3_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.2,1,.4,4]
d2 = [.05,.25,.1,1]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,2.8,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=4_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.2,1,0.5,5]
d2 = [.04,.2,.1,1]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,3.1,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=5_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.3,1.2,.6,6]
d2 = [.05,.2,.1,1]
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5)
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,3.1,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=6_0.5.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()



#r = 1.5
d1 = [.6,3,1.2,12]
d2 = [.4,2,.8,8]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.08,r'a)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=1.5_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


d1 = [.8,4,1.6,16]
d2 = [.4,2,.8,8]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.12,r'b)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=2_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.6,3.6,1.8,18]
d2 = [.2,1.2,.6,6]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',linewidth=3,label=ani[i], color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.21,r'c)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=3_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.8,4,1.6,16]
d2 = [.2,1,.4,4]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.3,r'd)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=4_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [.8,4,2,20]
d2 = [.16,.8,.4,4]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.32,r'e)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=5_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()

d1 = [1.2,4.8,2.4,24]
d2 = [.2,.8,.4,4]
for i in range(len(d1)):
    if i==1: continue
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle('../analysis/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+models[i]+'_'+mods[i]+'.pkl')
    H = data['ACF_h']
    V = data['ACF:_v']
    est=H/V
    s_H = data['sACF_h']
    s_V = data['sACF:_v']
    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    plt.plot(angles, est,'-',label=ani[i],linewidth=3, color = cols[i])
    plt.fill_between(angles, est-err, est+err,alpha=0.5, color = cols[i])
plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'autocorrelation ratio, $\alpha$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
plt.text(-10,1.33,r'f)',size=16)
plt.legend(prop={'size': 16}); 
plt.savefig('acf_systems_r=6_2.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]

plt.plot(R[0::4], fits[0::4],'s',color=cols[0],label='AP (fib.), '+r'$m=1$')
plt.plot(R[1::4], fits[1::4],'s',color=cols[1],label='AP (art.), '+r'$m=1$')
plt.plot(R[2::4], fits[2::4],'s',color=cols[2],label='GL, '+r'$m=1$')
plt.plot(R[3::4], fits[3::4],'s',color=cols[3],label='CB, '+r'$m=1$')

plt.plot(R[0::4], fits_2[0::4],'v',color=cols[0],label='AP (fib.), '+r'$m=1/2$')
plt.plot(R[1::4], fits_2[1::4],'v',color=cols[1],label='AP (art.), '+r'$m=1/2$')
plt.plot(R[2::4], fits_2[2::4],'v',color=cols[2],label='GL, '+r'$m=1/2$')
plt.plot(R[3::4], fits_2[3::4],'v',color=cols[3],label='CB, '+r'$m=1/2$')

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(Rm[0::3], fits_3[0::3],'^',color=cols[0],label='AP (fib.), '+r'$m=2$')
plt.plot(Rm[1::3], fits_3[1::3],'^',color=cols[2],label='GL, '+r'$m=2$')
plt.plot(Rm[2::3], fits_3[2::3],'^',color=cols[3],label='CB, '+r'$m=2$')

fits.extend(fits_2); fits.extend(fits_3) 
#R = np.array(R)
R.extend(R); R.extend(Rm)
c = np.polyfit(R,fits,1)
x = np.array([1,7])
plt.plot(x,x*c[0]+c[1],'k-',lw=3,label = 'Identity')

legend = axes.legend(fontsize=14,loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=3)

def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)
plt.xlim((1.3,6.2))
plt.grid()
plt.xlabel(r'r',size=16)
plt.ylabel(r'$r_{eff}$',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig('fittings.eps',format='eps',bbox_inches='tight')
plt.show()



fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]

plt.plot(R[0::4], np.array(fits_h[0::4])[:,0],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(R[1::4], np.array(fits_h[1::4])[:,0],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(R[2::4], np.array(fits_h[2::4])[:,0],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(R[3::4], np.array(fits_h[3::4])[:,0],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(R[0::4], np.array(fits_h_ii[0::4])[:,0],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[1::4], np.array(fits_h_ii[1::4])[:,0],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[2::4], np.array(fits_h_ii[2::4])[:,0],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(R[3::4], np.array(fits_h_ii[3::4])[:,0],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')


plt.plot(R[0::4], np.array(fits_v[0::4])[:,0],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(R[1::4], np.array(fits_v[1::4])[:,0],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(R[2::4], np.array(fits_v[2::4])[:,0],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(R[3::4], np.array(fits_v[3::4])[:,0],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(R[0::4], np.array(fits_v_ii[0::4])[:,0],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[1::4], np.array(fits_v_ii[1::4])[:,0],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[2::4], np.array(fits_v_ii[2::4])[:,0],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(R[3::4], np.array(fits_v_ii[3::4])[:,0],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(Rm[0::3], np.array(fits_h_iii[0::3])[:,0],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(Rm[1::3], np.array(fits_h_iii[1::3])[:,0],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(Rm[2::3], np.array(fits_h_iii[2::3])[:,0],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')

plt.plot(Rm[0::3], np.array(fits_v_iii[0::3])[:,0],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(Rm[1::3], np.array(fits_v_iii[1::3])[:,0],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(Rm[2::3], np.array(fits_v_iii[2::3])[:,0],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')


#fits.extend(fits_2); fits.extend(fits_3) 
#R = np.array(R)
'''R.extend(R); R.extend(Rm)
c = np.polyfit(R,fits,1)
x = np.array([1,7])
plt.plot(x,x*c[0]+c[1],'k-',lw=3,label = 'Linear fit')'''

#legend = axes.legend(fontsize=14,loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=3)

plt.xlim((1.3,6.2))
plt.grid()
plt.xlabel(r'r',size=16)
plt.ylabel(r'Multiplicative factor',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
'''def export_legend(legend, filename="legend.eps"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)'''
plt.text(.3,.7,r'a)',size=16)
plt.savefig('fit_params_a.eps',format='eps',bbox_inches='tight')
plt.show()



fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]

plt.plot(R[0::4], np.array(fits_h[0::4])[:,1],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(R[1::4], np.array(fits_h[1::4])[:,1],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(R[2::4], np.array(fits_h[2::4])[:,1],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(R[3::4], np.array(fits_h[3::4])[:,1],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(R[0::4], np.array(fits_h_ii[0::4])[:,1],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[1::4], np.array(fits_h_ii[1::4])[:,1],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[2::4], np.array(fits_h_ii[2::4])[:,1],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(R[3::4], np.array(fits_h_ii[3::4])[:,1],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')


plt.plot(R[0::4], np.array(fits_v[0::4])[:,1],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(R[1::4], np.array(fits_v[1::4])[:,1],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(R[2::4], np.array(fits_v[2::4])[:,1],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(R[3::4], np.array(fits_v[3::4])[:,1],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(R[0::4], np.array(fits_v_ii[0::4])[:,1],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[1::4], np.array(fits_v_ii[1::4])[:,1],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[2::4], np.array(fits_v_ii[2::4])[:,1],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(R[3::4], np.array(fits_v_ii[3::4])[:,1],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(Rm[0::3], np.array(fits_h_iii[0::3])[:,1],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(Rm[1::3], np.array(fits_h_iii[1::3])[:,1],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(Rm[2::3], np.array(fits_h_iii[2::3])[:,1],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')

plt.plot(Rm[0::3], np.array(fits_v_iii[0::3])[:,1],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(Rm[1::3], np.array(fits_v_iii[1::3])[:,1],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(Rm[2::3], np.array(fits_v_iii[2::3])[:,1],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')


#fits.extend(fits_2); fits.extend(fits_3) 
#R = np.array(R)
'''R.extend(R); R.extend(Rm)
c = np.polyfit(R,fits,1)
x = np.array([1,7])
plt.plot(x,x*c[0]+c[1],'k-',lw=3,label = 'Linear fit')'''

#legend = axes.legend(fontsize=14,loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=3)

plt.xlim((1.3,6.2))
plt.grid()
plt.xlabel(r'r',size=16)
plt.ylabel(r'Exponent (b,d)',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
'''def export_legend(legend, filename="legend.eps"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)'''
plt.text(.3,-.13,r'd)',size=16)
plt.savefig('fit_params_d_h.eps',format='eps',bbox_inches='tight')
plt.show()

fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]


plt.plot(np.array(fits_h[0::4])[:,1], np.array(fits_v[0::4])[:,1],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(np.array(fits_h[1::4])[:,1], np.array(fits_v[1::4])[:,1],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(np.array(fits_h[2::4])[:,1], np.array(fits_v[2::4])[:,1],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(np.array(fits_h[3::4])[:,1], np.array(fits_v[3::4])[:,1],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(np.array(fits_h_ii[0::4])[:,1], np.array(fits_v_ii[0::4])[:,1],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[1::4])[:,1], np.array(fits_v_ii[1::4])[:,1],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[2::4])[:,1], np.array(fits_v_ii[2::4])[:,1],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[3::4])[:,1], np.array(fits_v_ii[3::4])[:,1],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(np.array(fits_h_iii[0::3])[:,1], np.array(fits_v_iii[0::3])[:,1],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(np.array(fits_h_iii[1::3])[:,1], np.array(fits_v_iii[1::3])[:,1],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(np.array(fits_h_iii[2::3])[:,1], np.array(fits_v_iii[2::3])[:,1],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')
plt.plot([-0.22,-0.15],[-0.22,-0.15],'k-')
plt.xlim((-.225,-.15))
plt.grid()
plt.xlabel(r'b',size=16)
plt.ylabel(r'd',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.text(-.24,-.15,r'c)',size=16)
plt.savefig('fit_params_c_h.eps',format='eps',bbox_inches='tight')
plt.show()


fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]


plt.plot(np.array(fits_h[0::4])[:,0], np.array(fits_v[0::4])[:,0],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(np.array(fits_h[1::4])[:,0], np.array(fits_v[1::4])[:,0],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(np.array(fits_h[2::4])[:,0], np.array(fits_v[2::4])[:,0],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(np.array(fits_h[3::4])[:,0], np.array(fits_v[3::4])[:,0],'s',color=cols[3],label='CBM, '+r'$f=1$')

plt.plot(np.array(fits_h_ii[0::4])[:,0], np.array(fits_v_ii[0::4])[:,0],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[1::4])[:,0], np.array(fits_v_ii[1::4])[:,0],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[2::4])[:,0], np.array(fits_v_ii[2::4])[:,0],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(np.array(fits_h_ii[3::4])[:,0], np.array(fits_v_ii[3::4])[:,0],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(np.array(fits_h_iii[0::3])[:,0], np.array(fits_v_iii[0::3])[:,0],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(np.array(fits_h_iii[1::3])[:,0], np.array(fits_v_iii[1::3])[:,0],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(np.array(fits_h_iii[2::3])[:,0], np.array(fits_v_iii[2::3])[:,0],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')
plt.plot([0.3,0.7],[0.3,0.7],'k-')
#plt.xlim((1.3,6.2))
plt.grid()
plt.xlabel(r'a',size=14)
plt.ylabel(r'c',size=14)
plt.tick_params(axis='both', which='major', labelsize=14)

#plt.text(.21,.7,r'b)',size=16)
axes.text(-0.15, 1.05, r'b)', transform=axes.transAxes,
fontsize=14, fontweight='normal', va='top', ha='left')
plt.savefig('fit_params_b_h.eps',format='eps',bbox_inches='tight')
plt.show()


fig, axes = plt.subplots(1, 1)

R = [1.5, 1.5, 1.5, 1.5,
     2, 2, 2, 2, 3, 3, 3, 3,
     4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6]

R_to_save = []
b_to_save = []
model_to_save = []

#fits_rate = -np.array(fits_rate)
#fits_rate_ii = -np.array(fits_rate_ii)
#fits_rate_iii = -np.array(fits_rate_iii)

plt.plot(R[0::4], np.array(fits_rate[0::4])[:,0],'s',color=cols[0],label='AP (fib.), '+r'$f=1$')
plt.plot(R[1::4], np.array(fits_rate[1::4])[:,0],'s',color=cols[1],label='AP (art.), '+r'$f=1$')
plt.plot(R[2::4], np.array(fits_rate[2::4])[:,0],'s',color=cols[2],label='GLE, '+r'$f=1$')
plt.plot(R[3::4], np.array(fits_rate[3::4])[:,0],'s',color=cols[3],label='CBM, '+r'$f=1$')

R_to_save.extend(R[0::4])
b_to_save.extend(np.array(fits_rate[0::4])[:,0])
model_to_save.extend(['AP_chaos' for i in range(len(R[0::4]))])

R_to_save.extend(R[1::4])
b_to_save.extend(np.array(fits_rate[1::4])[:,0])
model_to_save.extend(['AP_ulrich' for i in range(len(R[1::4]))])

R_to_save.extend(R[2::4])
b_to_save.extend(np.array(fits_rate[2::4])[:,0])
model_to_save.extend(['GLE' for i in range(len(R[2::4]))])

R_to_save.extend(R[3::4])
b_to_save.extend(np.array(fits_rate[3::4])[:,0])
model_to_save.extend(['CBM' for i in range(len(R[3::4]))])

plt.plot(R[0::4], np.array(fits_rate_ii[0::4])[:,0],'v',color=cols[0],label='AP (fib.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[1::4], np.array(fits_rate_ii[1::4])[:,0],'v',color=cols[1],label='AP (art.), '+r'$f=1/\sqrt{2}$')
plt.plot(R[2::4], np.array(fits_rate_ii[2::4])[:,0],'v',color=cols[2],label='GLE, '+r'$f=1/\sqrt{2}$')
plt.plot(R[3::4], np.array(fits_rate_ii[3::4])[:,0],'v',color=cols[3],label='CBM, '+r'$f=1/\sqrt{2}$')

R_to_save.extend(R[0::4])
b_to_save.extend(np.array(fits_rate_ii[0::4])[:,0])
model_to_save.extend(['AP_chaos' for i in range(len(R[0::4]))])

R_to_save.extend(R[1::4])
b_to_save.extend(np.array(fits_rate_ii[1::4])[:,0])
model_to_save.extend(['AP_ulrich' for i in range(len(R[1::4]))])

R_to_save.extend(R[2::4])
b_to_save.extend(np.array(fits_rate_ii[2::4])[:,0])
model_to_save.extend(['GLE' for i in range(len(R[2::4]))])

R_to_save.extend(R[3::4])
b_to_save.extend(np.array(fits_rate_ii[3::4])[:,0])
model_to_save.extend(['CBM' for i in range(len(R[3::4]))])

Rm = [1.5, 1.5, 1.5,
     2, 2, 2, 3, 3, 3,
     4, 4, 4, 5, 5, 5,  6, 6, 6]

plt.plot(Rm[0::3], np.array(fits_rate_iii[0::3])[:,0],'^',color=cols[0],label='AP (fib.), '+r'$f=2\sqrt{2}$')
plt.plot(Rm[1::3], np.array(fits_rate_iii[1::3])[:,0],'^',color=cols[2],label='GLE, '+r'$f=2\sqrt{2}$')
plt.plot(Rm[2::3], np.array(fits_rate_iii[2::3])[:,0],'^',color=cols[3],label='CBM, '+r'$f=2\sqrt{2}$')


R_to_save.extend(Rm[0::3])
b_to_save.extend(np.array(fits_rate_iii[0::3])[:,0])
model_to_save.extend(['AP_chaos' for i in range(len(Rm[0::3]))])

R_to_save.extend(Rm[1::3])
b_to_save.extend(np.array(fits_rate_iii[1::3])[:,0])
model_to_save.extend(['GLE' for i in range(len(Rm[1::3]))])

R_to_save.extend(Rm[2::3])
b_to_save.extend(np.array(fits_rate_iii[2::3])[:,0])
model_to_save.extend(['CBM' for i in range(len(Rm[2::3]))])

data_to_save = np.c_[model_to_save,b_to_save, R_to_save]
data_to_save = np.array(data_to_save)
data_to_save = pd.DataFrame(data_to_save,columns=['model','b','r'])
data_to_save.to_csv('exponents_L='+str(L)+'_lag='+str(lag)+'.csv')

#fits_rate = list(fits_rate)
#fits_rate_ii = list(fits_rate_ii)
#fits_rate_iii = list(fits_rate_iii)

F = fits_rate
F.extend(fits_rate_ii); F.extend(fits_rate_iii); 
plt.plot([1.3,6.2],[np.mean(F),np.mean(F)],'k--',lw=3)

#fits.extend(fits_2); fits.extend(fits_3) 
#R = np.array(R)
'''R.extend(R); R.extend(Rm)
c = np.polyfit(R,fits,1)
x = np.array([1,7])
plt.plot(x,x*c[0]+c[1],'k-',lw=3,label = 'Linear fit')'''

#legend = axes.legend(fontsize=14,loc='upper center', bbox_to_anchor=(0.5, 1.5),ncol=3)

plt.xlim((1.3,6.2))
plt.grid()
plt.xlabel(r'r',size=14)
plt.ylabel(r'$b$',size=14)
plt.tick_params(axis='both', which='major', labelsize=14)
'''def export_legend(legend, filename="legend.eps"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)'''
#plt.text(.5,.26,r'b)',size=14)
axes.text(-0.2, 1.05, r'b)', transform=axes.transAxes,
fontsize=14, fontweight='normal', va='top', ha='left')
plt.savefig('fit_params_b.eps',format='eps',bbox_inches='tight')
plt.show()
