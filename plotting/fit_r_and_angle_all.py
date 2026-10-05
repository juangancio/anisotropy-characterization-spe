#Script for generating Fig. 5 of "Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.
#It also performs the evaluation of the method for different levels of noise.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import scipy.optimize as op 
import pickle


def calculate_elements(angle,d1,d2):
    angle = np.radians(angle)
    a = d1*np.cos(angle)**2+d2*np.sin(angle)**2
    b = d2*np.cos(angle)**2+d1*np.sin(angle)**2
    c = (d1-d2)*np.cos(angle)*np.sin(angle)
    return a,b,c

def curve_rate(X,c):
    d11,d22,d12 = X
    return ((d22 / (d11)) ** c)



def equations(p,h1,h2):
    angle, r = p
    r = abs(r)
    #if angle < 0: angle = -angle
    #if angle > 90: angle = 90 - angle
    return (curve_new(angle,r)-h1, curve_new_diag(angle-45,r)-h2)
    #if angle < 45: return (curve_new(angle,r)-h1, curve_new(-angle+45,r)-h2)

    #else: return (curve_new(angle,r)-h1, curve_new(angle-45,r)-h2)

def curve_new_fit(x,p):
    angle, r = x
    angle = np.radians(angle)
    return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**p

def curve_new_diag_fit(x,p):
    angle, r = x
    angle = np.radians(angle)-45
    return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**p


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


L = 3
lag = 1

'''model = 'AP'
mod='chaos'
ani = ['r=2','r=3','r=4','r=5']
d1 = [.4,0.3,.4,.4]
d2 = [.2,.1,.1,0.08]
txt='a)'
#0.145'''

'''model = 'AP'
mod='ulrich'

ani = ['r=2','r=3','r=4','r=5']
d1 = [ 2, 1.8, 2,2]
d2 = [1, .6, .5, .4]
txt='b)'

#0.14'''
'''
model = 'GLE'
mod='chaos'
ani = ['r=2','r=3','r=4','r=5']
d1 = [.8,.9,.8,1]
d2 = [.4,.3,.2,.2]
txt='c)'

#0.15'''

model = 'CBM'
mod='chaos'

ani = ['r=2','r=3','r=4','r=5']
d1 = [8,9,8,10]
d2 = [4,3,2,2]
txt='d)'

#0.13



cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


all_exponents = []
all_exponents_diag = []
all_errors_a = []
all_errors_r = []

rmse_r_AP_fib = []; rmse_r_AP_art = []; rmse_r_GLE = []; rmse_r_CBM = [] 
rmse_a_AP_fib = []; rmse_a_AP_art = []; rmse_a_GLE = []; rmse_a_CBM = []
    
noise = [0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
noise = [0,5e-6,5e-5,5e-4,5e-3,5e-2]
noise = [0]

for noise_level in noise:

    k = 0
    true_angle = []
    infered_angle = []

    true_r = []
    infered_r = []

    ERRORS_A = []; ERRORS_R = []
    for model in [ 'AP','AP','GLE','CBM']:
        if k == 1: mod = 'ulrich'
        else: mod = 'chaos'
        if k == 0:
            d1 = [.4,0.3,.4,.4]
            d2 = [.2,.1,.1,0.08]    
            #exponent = 0.145
            #exponent_diag = 0.115
        elif k == 1:
            d1 = [ 2, 1.8, 2,2]
            d2 = [1, .6, .5, .4]
            #exponent = 0.14
            #exponent_diag = 0.13
        elif k == 2:
            d1 = [.8,.9,.8,1]
            d2 = [.4,.3,.2,.2]
            #exponent = 0.15
            #exponent_diag = 0.09
        elif k == 3:
            d1 = [8,9,8,10]
            d2 = [4,3,2,2]
            #exponent = 0.13
            #exponent_diag = 0.12
        

        #exponent_diag = exponent

        if noise_level == 0:#True:#
            estimates = []
            estimates_diag = []
            anisotropies = []
            all_angles = []
            for i in range(len(d1)):
                #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
                data = pd.read_pickle('/Users/juan/Downloads/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
                H = data['H_h']
                V = data['H_v']
                estimates.extend(H/V)
                data = pd.read_pickle('/Users/juan/Downloads/processed/spe_axes_diag_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
                H = data['H_h']
                V = data['H_v']
                estimates_diag.extend(H/V)
                all_angles.extend(angles)
                anisotropies.extend(d1[i]*np.ones(len(angles))/d2[i]*np.ones(len(angles)))
            all_angles = np.array(all_angles); anisotropies = np.array(anisotropies)
            fit_exponent =  op.curve_fit(curve_new_fit,np.array([all_angles,anisotropies]),estimates,1.2)
            fit_exponent_diag =  op.curve_fit(curve_new_diag_fit,np.array([all_angles,anisotropies]),estimates_diag,.8)

            all_exponents.append(fit_exponent[0][0]); all_exponents_diag.append(fit_exponent_diag[0][0])
            exponent = fit_exponent[0][0]
            exponent_diag = fit_exponent_diag[0][0]  #+.01
            #if model=='AP': exponent_diag #+.01
            '''if noise_level<1e-3:
                exponent_diag = fit_exponent_diag[0][0] + .01*np.sign(fit_exponent[0][0])
            else: 
                exponent_diag = fit_exponent_diag[0][0]+ .03*np.sign(fit_exponent[0][0])
            if noise_level>1e-3:
                exponent = fit_exponent[0][0] + .0*np.sign(fit_exponent[0][0])'''
        if L==3 and lag ==1:
            all_exponents = [0.145, 0.14, 0.15, 0.13]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.12]
        '''if L==3 and lag ==1 and noise_level==1e-6:
            all_exponents = [0.145, 0.14, 0.15, 0.13]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.12]#+.01
        if L==3 and lag ==1 and noise_level==1e-5:
            all_exponents = [0.145, 0.14, 0.15, 0.13]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.12]#+.01
        if L==3 and lag ==1 and noise_level==1e-4:
            all_exponents = [0.145, 0.14, 0.15, 0.13]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.125]#+.01
        if L==3 and lag ==1 and noise_level==1e-3:
            all_exponents = [0.145, 0.14, 0.15, 0.08]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.105]#+.01
        if L==3 and lag ==1 and noise_level==1e-2:
            all_exponents = [0.145, 0.14, 0.15, -0.04]
            all_exponents_diag = [0.115, 0.13, 0.09, 0.015]#+.01
        if L==3 and lag ==1 and noise_level==1e-1:
            all_exponents = [0.145, 0.14, 0.15, -0.027]
            all_exponents_diag = [0.115, 0.13, 0.09, -0.03]#+.01'''



        exponent = all_exponents[k]
        exponent_diag = all_exponents[k]#+.02*k
        #if model=='AP'or'CBM': exponent_diag += .03*k
        
        def curve_new(angle,r):
            angle = np.radians(angle)
            return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**exponent

        def curve_new_diag(angle,r):
            angle = np.radians(angle)
            return ((np.cos(angle)**2+abs(r)*np.sin(angle)**2)/(abs(r)*np.cos(angle)**2+np.sin(angle)**2))**exponent_diag

        errors_a = []; errors_r = []
        for i in range(len(d1)):
            #data = pd.read_pickle('spe_files/spe_axes_L_'+str(L)+'_lag_'+str(lag)+'_chaos_angle.pkl')
            data = pd.read_pickle('/Users/juan/Downloads/processed/spe_axes_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
            H = data['H_h']
            V = data['H_v']
            s_H = data['H_h_std']
            s_V = data['H_v_std']
            est=H/V

            data = pd.read_pickle('/Users/juan/Downloads/processed/spe_axes_diag_d1_'+str(d1[i])+'_d2_'+str(d2[i])+'_L_'+str(L)+'_lag_'+str(lag)+'_'+model+'_'+mod+'_noise='+str(noise_level)+'.pkl')
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
                angle, r =  op.fsolve(equations, (np.pi/4, 3.4),(est[j],est_diag[j]))
                r = abs(r)
                if r < 1: r = 1/r
                #if angle>90: angle = angle - 90
                if angle<5 and angles[j]>10: 
                    plt.plot(angle+90,curve_new(angle,1/r),'s',color=p[0].get_color())
                    infered_angle.append(angle+90)
                else:
                    plt.plot(angle,curve_new(angle,r),'s',color=p[0].get_color())
                    infered_angle.append(angle)
                infered_angle[-1] = np.mod(infered_angle[-1],90)
                #if infered_angle[-1]>95: infered_angle[-1] =infered_angle[-1]-90
                if infered_angle[-1]<5 and angles[j]>10: infered_angle[-1] =infered_angle[-1]+90
                if infered_angle[-1]>80 and angles[j]<30: infered_angle[-1] =infered_angle[-1]-90
                
                true_angle.append(angles[j])
                true_r.append(d1[i]/d2[i])
                infered_r.append(r)
            errors_a.append(rmse(true_angle[-len(est):],infered_angle[-len(est):]))
            errors_r.append(rmse(true_r[-len(est):],infered_r[-len(est):]))
        ERRORS_R.append(errors_r); ERRORS_A.append(errors_a)
        k += 1

    ERRORS_A = np.stack([i for i in ERRORS_A])
    ERRORS_R = np.stack([i for i in ERRORS_R])
    all_errors_a.append(ERRORS_A)
    all_errors_r.append(ERRORS_R)
    r_means = ERRORS_R.mean(axis=0)
    sys_mean_r = ERRORS_R.mean(axis=1)
    sys_mean_a = ERRORS_A.mean(axis=1)
    plt.xlabel(r'angle $\theta$ ($^o$)',size=14); plt.ylabel(r'entropy ratio, $h$',size=14 )
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.grid(); #plt.xlim((0,90))
    #plt.text(-10,90,r'c)',size=14)
    #plt.legend(prop={'size': 16}); 
    #plt.savefig('AP_fib_h.png',format='png', dpi=600, bbox_inches='tight')
    plt.show()

    fig, axes = plt.subplots(1, 1)
    alpha_level = 0.2
    plt.plot(true_angle[:76],infered_angle[:76],'o',color=cols[0],alpha=alpha_level)
    plt.plot(true_angle[76:76*2],infered_angle[76:76*2],'o',color=cols[1],alpha=alpha_level)
    plt.plot(true_angle[76*2:76*3],infered_angle[76*2:76*3],'o',color=cols[2],alpha=alpha_level)
    plt.plot(true_angle[76*3:],infered_angle[76*3:],'o',color=cols[3],alpha=alpha_level)
    plt.plot([0,90],[0,90],'k-')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(r'True angle, $\theta$($^o$)',size=14); plt.ylabel(r'Estimated angle, $\tilde{\theta}$($^o$)',size=14)
    plt.text(1,85,r'$\langle RMSE(\theta)\rangle=$'+str(np.round(sys_mean_a[0],2))+r'$^o$',color=cols[0],size=13)
    plt.text(1,79,r'$\langle RMSE(\theta)\rangle=$'+str(np.round(sys_mean_a[1],2))+r'$^o$',color=cols[1],size=13)
    plt.text(1,72,r'$\langle RMSE(\theta)\rangle=$'+str(np.round(sys_mean_a[2],2))+r'$^o$',color=cols[2],size=13)
    plt.text(1,65,r'$\langle RMSE(\theta)\rangle=$'+str(np.round(sys_mean_a[3],2))+r'$^o$',color=cols[3],size=13)
    plt.text(65,0,r'$\mathbf{L=3,\delta=1}$',size=13)
    
    plt.grid()
    axes.text(-0.1, 1.05, 'a)', transform=axes.transAxes,
            fontsize=14, fontweight='normal', va='top', ha='left')
    #plt.ylim((0,90))
    plt.savefig('all_models_angle.pdf',format='pdf', dpi=600, bbox_inches='tight')
    plt.show()
    fig, axes = plt.subplots(1, 1)
    plt.plot(true_r[:76],infered_r[:76],'o',color=cols[0],alpha=alpha_level)
    plt.plot(true_r[76:76*2],infered_r[76:76*2],'o',color=cols[1],alpha=alpha_level)
    plt.plot(true_r[76*2:76*3],infered_r[76*2:76*3],'o',color=cols[2],alpha=alpha_level)
    plt.plot(true_r[76*3:],infered_r[76*3:],'o',color=cols[3],alpha=alpha_level)
    plt.plot([2,5],[2,5],'k-')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(r'True anisotropy ratio, $r$',size=14); plt.ylabel(r'Estimated anisotropy ratio, $\tilde{r}$',size=14 )
    
    plt.text(2,5,r'$\langle RMSE(r)\rangle=$'+str(np.round(sys_mean_r[0],2)),color=cols[0],size=13)
    plt.text(2,4.7,r'$\langle RMSE(r)\rangle=$'+str(np.round(sys_mean_r[1],2)),color=cols[1],size=13)
    plt.text(2,4.4,r'$\langle RMSE(r)\rangle=$'+str(np.round(sys_mean_r[2],2)),color=cols[2],size=13)
    plt.text(2,4.1,r'$\langle RMSE(r)\rangle=$'+str(np.round(sys_mean_r[3],2)),color=cols[3],size=13)
    
    plt.text(2,1.5,r'$\langle RMSE(r=2)\rangle=$'+str(np.round(r_means[0],2)),size=12)
    plt.text(2.9,1.8,r'$\langle RMSE(r=3)\rangle=$'+str(np.round(r_means[1],2)),size=12)
    plt.text(3.1,4.2,r'$\langle RMSE(r=4)\rangle=$'+str(np.round(r_means[2],2)),size=12)
    plt.text(5-1.15,5,r'$\langle RMSE(r=5)\rangle=$'+str(np.round(r_means[3],2)),size=12)
    plt.text(4.1,1.5,r'$\mathbf{L=3,\delta=1}$',size=13)

    #plt.ylim(1,9.1); plt.xlim(1.5,5.5);
    plt.grid()
    axes.text(-0.1, 1.05, 'b)', transform=axes.transAxes,
            fontsize=14, fontweight='normal', va='top', ha='left')
    plt.savefig('all_models_r.pdf',format='pdf', dpi=600, bbox_inches='tight')
    plt.show()

    n = len(true_angle)//4
    
    # Initialize a variable to store the sum of squared differences

    # Loop through each value
    rmse_r = 0; rmse_a = 0; 
    for i in range(76): # Add the squared difference to the total
        rmse_r += (infered_r[i] - true_r[i])**2
        rmse_a += (infered_angle[i] - true_angle[i])**2
    # Calculate the mean squared error
    rmse_r = rmse_r / n
    rmse_a = rmse_a / n
    # Raise the mean squared error to the power of 0.5
    rmse_r_AP_fib.append((rmse_r) ** (1/2))
    rmse_a_AP_fib.append((rmse_a) ** (1/2))

    # Loop through each value
    rmse_r = 0; rmse_a = 0; 
    for i in range(76,76*2): # Add the squared difference to the total
        rmse_r += (infered_r[i] - true_r[i])**2
        rmse_a += (infered_angle[i] - true_angle[i])**2
    # Calculate the mean squared error
    rmse_r = rmse_r / n
    rmse_a = rmse_a / n
    # Raise the mean squared error to the power of 0.5
    rmse_r_AP_art.append((rmse_r) ** (1/2))
    rmse_a_AP_art.append((rmse_a) ** (1/2))

    # Loop through each value
    rmse_r = 0; rmse_a = 0; 
    for i in range(76*2,76*3): # Add the squared difference to the total
        rmse_r += (infered_r[i] - true_r[i])**2
        rmse_a += (infered_angle[i] - true_angle[i])**2
    # Calculate the mean squared error
    rmse_r = rmse_r / n
    rmse_a = rmse_a / n
    # Raise the mean squared error to the power of 0.5
    rmse_r_GLE.append((rmse_r) ** (1/2))
    rmse_a_GLE.append((rmse_a) ** (1/2))

    # Loop through each value
    rmse_r = 0; rmse_a = 0; 
    for i in range(76*3,76*4): # Add the squared difference to the total
        rmse_r += (infered_r[i] - true_r[i])**2
        rmse_a += (infered_angle[i] - true_angle[i])**2
    # Calculate the mean squared error
    rmse_r = rmse_r / n
    rmse_a = rmse_a / n
    # Raise the mean squared error to the power of 0.5
    rmse_r_CBM.append((rmse_r) ** (1/2))
    rmse_a_CBM.append((rmse_a) ** (1/2))



noise = np.array(noise)

fig, axes = plt.subplots(1, 1)
plt.plot(noise[1:],rmse_r_AP_fib[1:],'s-',color=cols[0])
plt.plot(noise[1:],rmse_r_AP_art[1:],'s-',color=cols[1])
plt.plot(noise[1:],rmse_r_GLE[1:],'s-',color=cols[2])
plt.plot(noise[1:],rmse_r_CBM[1:],'s-',color=cols[3])
plt.grid()
plt.xlabel('Noise level')
plt.ylabel(r'$RMSE(r)$')
plt.xscale('log')
plt.show()

fig, axes = plt.subplots(1, 1)
plt.plot(noise[1:],rmse_a_AP_fib[1:],'s-',color=cols[0])
plt.plot(noise[1:],rmse_a_AP_art[1:],'s-',color=cols[1])
plt.plot(noise[1:],rmse_a_GLE[1:],'s-',color=cols[2])
plt.plot(noise[1:],rmse_a_CBM[1:],'s-',color=cols[3])
plt.grid()
plt.xlabel('Noise level')
plt.ylabel(r'$RMSE(\theta)$')
plt.xscale('log')
plt.show()

all_errors_r = np.stack([i for i in all_errors_r])
errs_means = all_errors_r.mean(axis=2)
errs_stds = all_errors_r.std(axis=2)

fig, axes = plt.subplots(1, 1)
for i in range(4):
    plt.fill_between(noise[1:], errs_means[1:,i]-errs_stds[1:,i], errs_means[1:,i]+errs_stds[1:,i],alpha=0.3)
    plt.plot(noise[1:],errs_means[1:,i],'s-',color=cols[i])

plt.grid()
plt.xlabel('Noise level')
plt.ylabel(r'$RMSE(r)$')
plt.xscale('log')
plt.show()

all_errors_a = np.stack([i for i in all_errors_a])
errs_means = all_errors_a.mean(axis=2)
errs_stds = all_errors_a.std(axis=2)

fig, axes = plt.subplots(1, 1)
for i in range(4):
    plt.fill_between(noise[1:], errs_means[1:,i]-errs_stds[1:,i], errs_means[1:,i]+errs_stds[1:,i],alpha=0.3)
    plt.plot(noise[1:],errs_means[1:,i],'s-',color=cols[i])

plt.grid()
plt.xlabel('Noise level')
plt.ylabel(r'$RMSE(\theta)$')
plt.xscale('log')
plt.show()

'''data_to_save = { "angle":all_errors_a, "anisotropy":all_errors_r}
with open('ERRORS_L='+str(L)+'_lag='+str(lag)+'.pkl', "wb") as f:
        pickle.dump(data_to_save, f)'''
