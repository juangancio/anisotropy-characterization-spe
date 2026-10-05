# Script for generating Fig. 7a of "Characterizing anisotropy in excitable media using spatial permutation entropy" (2026) by Gancio et al.
# It also plots Fig. 7 of the supplementary material.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as op 




def calculate_elements(df):
    a = df['d1']*np.cos(df['angle'])**2+df['d2']*np.sin(df['angle'])**2
    b = df['d2']*np.cos(df['angle'])**2+df['d1']*np.sin(df['angle'])**2
    c = (df['d1']-df['d2'])*np.cos(df['angle'])*np.sin(df['angle'])
    return a,b,c

def curve(X,c):
    d11,d22,d12 = X
    return ((d22/d11) ** c)
    #return (b*np.log(d11/d22)+c)
def curve_hor(X,b,c):
    d11,d22,d12 = X
    return (b*(  (d11)) ** c)
    #return (b*np.log(1/d11)+c)

def curve_ver(X,b,c):
    d11,d22,d12 = X
    return (b*(  (d22)) ** c)

def derivative(d,delta):
    return 0.213*d**(0.213-1)*delta

#L = 3; b_mean = 0.139
L = 4; b_mean = 0.213
#L = 5; b_mean = 0.262
raw_data = pd.read_csv('/Users/juan/downloads/rguimera-machine-scientist-c790f7cf35b5/test_L='+str(L)+'.csv')

filt_data = raw_data[raw_data['model'] != 'AP_chaos']
#filt_data = raw_data
#filt_data = filt_data[filt_data['d1']/filt_data['d2']  == 2]
x, y = filt_data['angle'], filt_data['H']/filt_data['V']

D11,D22,D12 = calculate_elements(filt_data)

fit = op.curve_fit(curve,(D11,D22,D12),y)
print(fit[0])


pred = curve((D11,D22,D12),fit[0][0])
plt.plot(y,D22/D11,'ro',label='r=4')
#plt.plot([.3,1],[.3,1],'k-')
plt.legend();plt.show()

r = [2,3,4,5]
cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
i=0
error = []
error_2 = []
fig, axes = plt.subplots(1, 1)
for model in ['AP_chaos','AP_ulrich','CBM_chaos','GLE_chaos']:
    param_1 = []; param_2 = []

    filt_data = raw_data[raw_data['model'] == model]
    D11,D22,D12 = calculate_elements(filt_data)

    '''x, y = filt_data['angle'], filt_data['H']/filt_data['V']
    fit = op.curve_fit(curve,(D11,D22,D12),y)
    pred = curve((D11,D22,D12),fit[0][0])
    error.append(sum((pred-y)**2))
    pred_2 = curve((D11,D22,D12),0.2)
    error_2.append(sum((pred_2-y)**2))'''

    

    x, y, z, f = filt_data['angle'], filt_data['H'], filt_data['V'], filt_data['mult']
    x, z = filt_data['angle'], filt_data['V']
    pred = curve_hor((D11,D22,D12),filt_data['fit1_hor'],filt_data['fit2_hor'])
    pred_z = curve_ver((D11,D22,D12),filt_data['fit1_ver'],filt_data['fit2_ver'])
    f = np.array(f); y = np.array(y); z = np.array(z)
    pred = np.array(pred);  pred_z = np.array(pred_z); 
    for j in range(len(f)):
        
        if f[j]<0.9: sym='v'
        elif f[j]>1.5: sym='^'
        else: sym='s'
        plt.plot(y[j],pred[j],sym,color=cols[i],alpha=.25)#,mfc='None')
        
        plt.plot(z[j],pred_z[j],sym,color=cols[i],alpha=.25)#,mfc='None')

    i+=1

#plt.plot([.7,1.45],[.7,1.45],'k-',lw=2)
plt.plot([.3,1.],[.3,1.],'k-')
plt.grid()
plt.ylabel(r'Approximation',size=16)
plt.xlabel(r'True',size=16)
plt.tick_params(axis='both', which='major', labelsize=16)
#plt.text(.5,1.6,r'a)',size=16)
# plt.text(.15,1,r'a)',size=16)

axes.text(-0.2, 1.01, 'a)', transform=axes.transAxes,
        fontsize=16, fontweight='normal', va='top', ha='left')

plt.savefig('fit_params_a_h.pdf',format='pdf',bbox_inches='tight')
plt.show()
plt.show()



i=0
error = []
error_2 = []
min_h = 1
max_h = 1
fig, axes = plt.subplots(1, 1)
for model in ['AP_chaos','AP_ulrich','CBM_chaos','GLE_chaos']:
    param_1 = []; param_2 = []

    filt_data = raw_data[raw_data['model'] == model]
    D11,D22,D12 = calculate_elements(filt_data)

    '''x, y = filt_data['angle'], filt_data['H']/filt_data['V']
    fit = op.curve_fit(curve,(D11,D22,D12),y)
    pred = curve((D11,D22,D12),fit[0][0])
    error.append(sum((pred-y)**2))
    pred_2 = curve((D11,D22,D12),0.2)
    error_2.append(sum((pred_2-y)**2))'''

    

    x, y, f = filt_data['angle'], filt_data['H']/filt_data['V'], filt_data['mult']
    fit = op.curve_fit(curve,(D11,D22,D12),y)
    pred = curve((D11,D22,D12),fit[0][0])
    error.append(sum((pred-y)**2))
    pred_2 = curve((D11,D22,D12),0.2)
    error_2.append(sum((pred_2-y)**2))
    f = np.array(f); y = np.array(y); 
    pred = np.array(pred);  pred_z = np.array(pred_z); 
    for j in range(len(f)):
        
        if f[j]<0.9: sym='v'
        elif f[j]>1.5: sym='^'
        else: sym='s'
        plt.plot(y[j],pred[j],sym,color=cols[i],alpha=.3)#,mfc='None')
        if min(y)<min_h:min_h=min(y)
        if max(y)>max_h:max_h=max(y)


    i+=1

plt.plot([.7,1.45],[.7,1.45],'k-',lw=2)
#plt.plot([.3,1.],[.3,1.],'k-')
plt.grid()
plt.ylabel(r'$\tilde{h}$',size=14)
plt.xlabel(r'$h$',size=14)
plt.tick_params(axis='both', which='major', labelsize=14)
#plt.text(.5,1.6,r'a)',size=16)
# plt.text(.15,1,r'a)',size=16)

axes.text(-0.2, 1.01, 'a)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
plt.text(1.3,.7,r'$\mathbf{L=4, \delta=1}$',size=16)
#plt.savefig('fit_params_a.pdf',format='pdf',bbox_inches='tight')
plt.show()
plt.show()


'''D=[]

for delta in np.arange(.01,.2,.01):
    diff = []
    for d in np.arange(0.5,6,.1):
        diff.append(abs(derivative(d,delta)))
    D.append(diff)
D = np.array(D)
plt.imshow(D,norm='log',aspect='auto',origin='lower',cmap='plasma',interpolation='bilinear')

plt.colorbar()'''


r = [2,3,4,5]
cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
i=0
E_mean = []
E_std = []
E2_mean = []
E2_std = []



for model in ['AP_chaos','AP_ulrich','CBM_chaos','GLE_chaos']:
    error = []
    error_2 = []
    for mult in [0.5,1,2]:
        if mult == 2 and model == 'AP_ulrich': continue
        for r in [1.5,2,3,4,5,6]:
            param_1 = []; param_2 = []

            filt_data = raw_data[raw_data['model'] == model]
            filt_data = filt_data[np.round(filt_data['r'],2) == r]
            filt_data = filt_data[filt_data['mult'] == mult]
            D11,D22,D12 = calculate_elements(filt_data)

            x, y = filt_data['angle'], filt_data['H']/filt_data['V']
            fit = op.curve_fit(curve,(D11,D22,D12),y)
            pred = curve((D11,D22,D12),fit[0][0])
            error.append(np.sqrt(sum((pred-y)**2)/19))
            pred_2 = curve((D11,D22,D12),b_mean)
            error_2.append(np.sqrt(sum((pred_2-y)**2)/19))
            #x, y = filt_data['angle'], filt_data['H']
            #pred = curve_hor((D11,D22,D12),filt_data['fit1_hor'],filt_data['fit2_hor'])

        
            #x, y = filt_data['angle'], filt_data['V']
            #pred = curve_ver((D11,D22,D12),filt_data['fit1_ver'],filt_data['fit2_ver'])
            #plt.plot(pred,y,'s',color=cols[i],alpha=.6,mfc='None')
            i+=1
    E_mean.append(np.mean(error)); E_std.append(np.std(error)/np.sqrt(len(error))); 
    E2_mean.append(np.mean(error_2)); E2_std.append(np.std(error_2)/np.sqrt(len(error_2))); 

barWidth = 0.3
r1 = np.arange(len(E_mean))
r2 = [x + barWidth for x in r1]
#fig, axes = plt.subplots(1, 1)
fig, axes = plt.subplots(1, 1, figsize=(8,3.6))
# Create blue bars
plt.grid(which='both',axis='y')
plt.bar(r1, E_mean, width = barWidth, color = cols[0], yerr=E_std,edgecolor = 'black', capsize=7, label='Fitted param.')
# Create cyan bars
plt.bar(r2, E2_mean, width = barWidth, color = cols[1],yerr=E2_std, edgecolor = 'black', capsize=7, label='Fixed param.')

# general layout
plt.xticks([r + barWidth for r in range(len(E_mean))], ['AP (fib)', 'AP (art.)', 'CB', 'GL'])
plt.ylabel(r'$\langle RMSE(h)\rangle$',fontsize=14); plt.xlabel('Model',fontsize=14)
plt.tick_params(axis='both', labelsize=14)
#plt.legend(fontsize=12, bbox_to_anchor=(0.2, 1),ncol=2)
#axes.text(-0.1, 1.05, 'c)', transform=axes.transAxes,fontsize=14, fontweight='normal', va='top', ha='left')
# Show graphic
axes.set_yscale('log')
axes.set_axisbelow(True)
axes.text(-0.13, 1.05, 'c)', transform=axes.transAxes,
        fontsize=14, fontweight='normal', va='top', ha='left')
axes.text(0.3,.85,r'$L=$'+str(L)+r'$, \delta=1$', transform=axes.transAxes,size=16)
#plt.ylim(0,.155)
plt.savefig('errors_L='+str(L)+'.pdf',format='pdf',bbox_inches='tight')
plt.show()
