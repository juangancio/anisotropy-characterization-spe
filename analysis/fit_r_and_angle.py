#Main script for implementing the algorithm proposed in
#"Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.,
#for simoulstanously estimating the orientation and anisotropy ratio.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import scipy.optimize as op 


def calculate_elements(angle,d1,d2):
    angle = np.radians(angle)
    a = d1*np.cos(angle)**2+d2*np.sin(angle)**2
    b = d2*np.cos(angle)**2+d1*np.sin(angle)**2
    c = (d1-d2)*np.cos(angle)*np.sin(angle)
    return a,b,c

def curve_rate(X,c):
    d11,d22,d12 = X
    return ((d22 / (d11)) ** c)

def curve_new(angle,r):
    angle = np.radians(angle)
    return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**0.15

def curve_new_diag(angle,r):
    angle = np.radians(angle)
    return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**0.125


def equations(p,h1,h2):
    angle, r = p

    #if angle < 0: angle = -angle
    #if angle > 90: angle = 90 - angle
    return (curve_new(angle,r)-h1, curve_new_diag(angle-45,r)-h2)
    #if angle < 45: return (curve_new(angle,r)-h1, curve_new(-angle+45,r)-h2)

    #else: return (curve_new(angle,r)-h1, curve_new(angle-45,r)-h2)

def rmse(true,infered):
    n = len(true)
    # Loop through each value
    rmse_r = 0; 
    for i in range(n): # Add the squared difference to the total
        rmse_r += (infered[i] - true[i])**2
    # Calculate the mean squared error
    rmse_r = rmse_r / n
    # Raise the mean squared error to the power of 0.5
    return (rmse_r) ** (1/2)



angles=np.linspace(0,90,19)


L=3
lag=1

'''model = 'AP'
mod='chaos'
ani = ['r=2','r=3','r=4','r=5']
d1 = [.4,0.3,.4,.4]
d2 = [.2,.1,.1,0.08]
txt='a)'
#0.145/0.115'''

'''model = 'AP'
mod='ulrich'

ani = ['r=2','r=3','r=4','r=5']
d1 = [ 2, 1.8, 2,2]
d2 = [1, .6, .5, .4]
txt='b)'''

#0.14/0.13

'''model = 'GLE'
mod='chaos'
ani = ['r=2','r=3','r=4','r=5']
d1 = [.8,.9,.8,1]
d2 = [.4,.3,.2,.2]
txt='c)'

#0.15/0.09'''

'''model = 'CBM'
mod='chaos'

ani = ['r=2','r=3','r=4','r=5']
d1 = [8,9,8,10]
d2 = [4,3,2,2]
txt='d)'''

#0.13/0.12

model = 'AP'
mod='spiral'

ani = ['r=2','r=3','r=4','r=5']
d1 = [ 1.4, 1.8, 2,2]
d2 = [.7, .6, .5, .4]
txt='b)'
#0.15/0.125


noise_level = 0

cols=[]

d1_to_save = []
d2_to_save = []
angles_to_save = []
entropy_to_save = []
r_to_save = []

path = '/Users/juan/Downloads/'
path = ''

fig, axes = plt.subplots(1, 1)

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle(path+'processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
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
    c = []
    for j in range(len(angles)): c.append(curve_new(angles[j],d1[i]/d2[i]))
    plt.plot(angles,c,'--',color='k')
    #if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    #else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    cols.append(p[0].get_color())
    '''d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)'''


plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); plt.xlim((0,90))
#plt.text(-10,1.3,r'a)',size=16)
axes.text(-0.15, 1.05, r'b)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.legend(prop={'size': 16}); 
plt.savefig(model+'_'+mod+'_h_vs_angle.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()
fig, axes = plt.subplots(1, 1)
for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle(path+'processed/spe_axes_diag_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V

    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p=plt.plot(angles-45, est,label=ani[i],lw=3)
    plt.fill_between(angles-45, est-err, est+err,alpha=0.5)
    type(p) # list
    type(p[0]) # <class 'matplotlib.lines.Line2D'>
    D11,D22,D12 = calculate_elements(angles-45,d1[i]*np.ones(angles.shape),d2[i]*np.ones(angles.shape))
    fit = op.curve_fit(curve_rate,(D11,D22,D12),est)
    c = []
    for j in range(len(angles)): c.append(curve_new_diag(angles[j]-45,d1[i]/d2[i]))
    plt.plot(angles-45,c,'--',color='k')
    #if i == len(d1)-1:plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k',label='fits')
    #else: plt.plot(angles,curve_rate((D11,D22,D12),fit[0][0]),'--',color='k')
    cols.append(p[0].get_color())
    '''d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)'''


plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); 
#plt.text(-20-45,1.02,r'b)',size=16)
axes.text(-0.15, 1.05, r'b)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
#plt.legend(prop={'size': 16}); 
plt.savefig(model+'_'+mod+'_h_vs_angle_diag.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

true_angle = []
infered_angle = []

true_r = []
infered_r = []
errors_a = []
errors_r = []

for i in range(len(d1)):
    #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
    data = pd.read_pickle(path+'processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est=H/V

    data = pd.read_pickle(path+'processed/spe_axes_diag_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
    H = data['H_h']
    V = data['H_v']
    s_H = data['H_h_std']
    s_V = data['H_v_std']
    est_diag=H/V

    err=np.sqrt((s_H/V)**2+(H/(V**2)*s_V)**2)
    p=plt.plot(angles, est,label=ani[i],lw=2)
    plt.fill_between(angles, est-err, est+err,alpha=0.3)
    type(p) # list
    type(p[0])

    for j in range(len(est)):
        angle, r =  op.fsolve(equations, (6*np.pi/8, 3.4),(est[j],est_diag[j]))
        r = abs(r)
        if r < 1: r = 1/r
        if angle>90: angle = angle - 90
        if angle<-5: angle = angle + 90
        if angle<5 and angles[j]>10: 
            plt.plot(angle+90,curve_new(angle,1/r),'s',color=p[0].get_color())
            infered_angle.append(angle+90)
        else:
            plt.plot(angle,curve_new(angle,r),'s',color=p[0].get_color())
            infered_angle.append(angle)
        if infered_angle[-1]>90: infered_angle[-1] =infered_angle[-1]-90
        if infered_angle[-1]<5 and angles[j]>10: infered_angle[-1] =infered_angle[-1]+90
        #if angle < -5 or angle >95: infered_angle[-1] =np.mod(infered_angle[-1],90)
      
        true_angle.append(angles[j])
        true_r.append(d1[i]/d2[i])
        infered_r.append(r)



     # <class 'matplotlib.lines.Line2D'>
    '''d1_to_save.extend(d1[i]*np.ones(angles.shape))
    d2_to_save.extend(d2[i]*np.ones(angles.shape))
    angles_to_save.extend(np.radians(angles))
    entropy_to_save.extend(est)'''

    errors_a.append(rmse(true_angle[-len(est):],infered_angle[-len(est):]))
    errors_r.append(rmse(true_r[-len(est):],infered_r[-len(est):]))

plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); #plt.xlim((0,90))
#plt.text(-10,90,r'c)',size=16)
#plt.legend(prop={'size': 16}); 
#plt.savefig('AP_fib_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()
fig, axes = plt.subplots(1, 1)
alpha_level = 0.5
plt.plot(true_angle[:19],infered_angle[:19],'s',color=cols[0],alpha=alpha_level)
plt.plot(true_angle[19:19*2],infered_angle[19:19*2],'s',color=cols[1],alpha=alpha_level)
plt.plot(true_angle[19*2:19*3],infered_angle[19*2:19*3],'s',color=cols[2],alpha=alpha_level)
plt.plot(true_angle[19*3:],infered_angle[19*3:],'s',color=cols[3],alpha=alpha_level)
plt.plot([0,90],[0,90],'k-')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r'true angle $\theta$ ($^o$)',size=16); plt.ylabel(r'estimated angle $\theta$ ($^o$)',size=16 )
plt.grid()
axes.text(-0.15, 1.05, 'c)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
#plt.ylim((0,90))
plt.savefig(model+'_'+mod+'_angle.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()
fig, axes = plt.subplots(1, 1)
plt.plot(true_r[:19],infered_r[:19],'s',color=cols[0],alpha=alpha_level)
plt.plot(true_r[19:19*2],infered_r[19:19*2],'s',color=cols[1],alpha=alpha_level)
plt.plot(true_r[19*2:19*3],infered_r[19*2:19*3],'s',color=cols[2],alpha=alpha_level)
plt.plot(true_r[19*3:],infered_r[19*3:],'s',color=cols[3],alpha=alpha_level)
plt.plot([2,5],[2,5],'k-')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r'true anisotropy ratio, $r$',size=16); plt.ylabel(r'estimated anisotropy ratio, $r$',size=16 )
plt.grid()
axes.text(-0.15, 1.05, 'd)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.savefig(model+'_'+mod+'_r.pdf',format='pdf', dpi=600, bbox_inches='tight')
plt.show()


plt.plot(true_angle[:19],infered_r[:19],'s',color=cols[0],alpha=alpha_level)
plt.plot(true_angle[19:19*2],infered_r[19:19*2],'s',color=cols[1],alpha=alpha_level)
plt.plot(true_angle[19*2:19*3],infered_r[19*2:19*3],'s',color=cols[2],alpha=alpha_level)
plt.plot(true_angle[19*3:],infered_r[19*3:],'s',color=cols[3],alpha=alpha_level)
plt.plot([2,5],[2,5],'k-')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r'true anisotropy ratio, $r$',size=16); plt.ylabel(r'estimated anisotropy ratio, $r$',size=16 )
plt.grid()
plt.text(1.5,5.25,r'd)',size=16)
#plt.savefig(model+'_'+mod+'_r.png',format='png', dpi=600, bbox_inches='tight')
plt.show()

'''
angles=np.linspace(0,90,19)-45

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
    data = pd.read_pickle('processed/spe_axes_diag_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'.pkl')
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
    #d1_to_save.extend(d1[i]*np.ones(angles.shape))
    #d2_to_save.extend(d2[i]*np.ones(angles.shape))
    #angles_to_save.extend(np.radians(angles))
    #entropy_to_save.extend(est


plt.xlabel(r'angle $\theta$ ($^o$)',size=16); plt.ylabel(r'entropy ratio, $h$',size=16 )
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(); #plt.xlim((0,90))
#plt.text(-10,1.55,r'a)',size=16)
plt.legend(prop={'size': 16}); 
#plt.savefig('AP_fib_h.png',format='png', dpi=600, bbox_inches='tight')
plt.show()'''
i=17
x, y =  op.fsolve(equations, (np.pi/4, 3),(est[i],est_diag[i]))
